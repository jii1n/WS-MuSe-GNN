import numpy as np
import torch
import os
from torch_geometric.data import DataLoader, Data
from tqdm import tqdm
from split import StratifiedGroupKFold
from pathlib import Path
from sklearn.model_selection import StratifiedShuffleSplit
import pandas as pd
from data.gene_expression_dataset import DatasetFromCSV
from data.gene_interaction_graph import StringDBGraph
from data.gene_expression_graph import GeneExpressionGraph
from sklearn.metrics import f1_score,classification_report



def cross_validate(model, X, labels, device, config, experiments, experiment_name,
                   num_folds, learning_rate, batch_size, weight_decay, epochs,
                   mixsplit=False):
    Path(f"gnn_checkpoints").mkdir(exist_ok=True)
    if mixsplit:
        stratified_k_fold = StratifiedShuffleSplit(num_folds, test_size=0.1, random_state=config.seed)
    else:
        stratified_k_fold = StratifiedGroupKFold(num_folds)

    device = device
    metrics = {
        'all_val_accs': [],
        'all_losses': [],
        'all_f1_scores': {'micro': [], 'macro': [], 'weighted': []},
        'class_metrics': {'long-lived': [], 'normal-lived': [], 'short-lived': []}
    }
    fold_predictions = []
    val_indices_list = []
    val_folds_data = [] 

    for fold, (train_index, val_index) in enumerate(stratified_k_fold.split(np.zeros(len(labels)), labels.values.reshape(-1), experiments)):
        val_indices_list.append(val_index)
        Path(f"/data/bi1/results/neural/{experiment_name}").mkdir(parents=True, exist_ok=True)
        f = open(f"/data/bi1/results/neural/{experiment_name}/fold_{fold}_training_progress.txt", "w")
        model.apply(model.weight_reset)

        if config.train_GNN:
            xy_train_fold = X[list(train_index)]
            xy_val_fold = X[list(val_index)]
        else:
            xy_train_fold = [Data(x=torch.Tensor(x.values), y=labels.loc[idx].values) for idx, x in X.iloc[list(train_index)].iterrows()]
            xy_val_fold = [Data(x=torch.Tensor(x.values), y=labels.loc[idx].values) for idx, x in X.iloc[list(val_index)].iterrows()]
        val_folds_data.append(xy_val_fold)

        opt = torch.optim.Adam(model.parameters(), lr=learning_rate, weight_decay=weight_decay)
        train_loader = DataLoader(xy_train_fold, batch_size=batch_size, shuffle=True)
        val_loader = DataLoader(xy_val_fold, batch_size=batch_size, shuffle=True)

        fold_val_accs = []
        fold_losses = []
        fold_f1_scores_micro = []
        fold_f1_scores_macro = []
        fold_f1_scores_weighted = []
        all_epoch_predictions = []
        loaded_upper_checkpoint = False

        for epoch in range(1, epochs + 1):
            upper_checkpoint = int((epoch + config.eval_model_every) / config.eval_model_every) * config.eval_model_every if epoch % config.eval_model_every != 0 and 
epoch != 1 else epoch
            upper_checkpoint_path = f"gnn_checkpoints/fold_{fold}_epoch_{upper_checkpoint}_{experiment_name}.pt.tar"
            path = f"gnn_checkpoints/fold_{fold}_epoch_{epoch}_{experiment_name}.pt.tar"

            if os.path.exists(upper_checkpoint_path):
                if not loaded_upper_checkpoint:
                    checkpoint = torch.load(upper_checkpoint_path)
                    model.load_state_dict(checkpoint['model_state_dict'])
                    opt.load_state_dict(checkpoint['optimizer_state_dict'])
                    total_loss = checkpoint['loss']
                    loaded_upper_checkpoint = True
            else:
                model, opt, total_loss = train_epoch(train_loader, model, opt, device)

            if epoch % config.eval_model_every == 0 or epoch == 1:
                if config.train_GNN and not loaded_upper_checkpoint:
                    torch.save({
                        'epoch': epoch,
                        'model_state_dict': model.state_dict(),
                        'optimizer_state_dict': opt.state_dict(),
                        'loss': total_loss,
                    }, path)
                loaded_upper_checkpoint = False

                val_acc, epoch_predictions, f1_scores, class_f1_scores = test_epoch(val_loader, model, device)
                y_true = [data.y.item() for data in xy_val_fold]
                y_pred = epoch_predictions

                for class_name, class_index in zip(['long-lived', 'normal-lived', 'short-lived'], [0, 1, 2]):
                    class_indices = [i for i, val in enumerate(y_true) if val == class_index]
                    class_y_true = [y_true[i] for i in class_indices]
                    class_y_pred = [y_pred[i] for i in class_indices]
                    class_f1 = f1_score(class_y_true, class_y_pred, average=None, zero_division=0)
                    metrics['class_metrics'][class_name].append(np.mean(class_f1))

                overall_micro_f1 = f1_score(y_true, y_pred, average='micro', zero_division=0)
                overall_macro_f1 = f1_score(y_true, y_pred, average='macro', zero_division=0)
                overall_weighted_f1 = f1_score(y_true, y_pred, average='weighted', zero_division=0)

                fold_val_accs.append(val_acc)
                fold_losses.append(total_loss)
                fold_f1_scores_micro.append(overall_micro_f1)
                fold_f1_scores_macro.append(overall_macro_f1)
                fold_f1_scores_weighted.append(overall_weighted_f1)
                all_epoch_predictions.append(epoch_predictions)

                f.write(f"Fold {fold}. Epoch {epoch}. Loss: {total_loss}. Validation accuracy: {val_acc}. \n")

        metrics['all_val_accs'].append(fold_val_accs)
        metrics['all_losses'].append(fold_losses)
        metrics['all_f1_scores']['micro'].append(fold_f1_scores_micro)
        metrics['all_f1_scores']['macro'].append(fold_f1_scores_macro)
        metrics['all_f1_scores']['weighted'].append(fold_f1_scores_weighted)
        fold_predictions.append(all_epoch_predictions)

    avg_val_accs_across_folds = [float(sum(col)) / len(col) for col in zip(*metrics['all_val_accs'])]
    avg_loss_across_folds = [float(sum(col)) / len(col) for col in zip(*metrics['all_losses'])]
    avg_f1_scores_micro = [float(sum(col)) / len(col) for col in zip(*metrics['all_f1_scores']['micro'])]
    avg_f1_scores_macro = [float(sum(col)) / len(col) for col in zip(*metrics['all_f1_scores']['macro'])]
    avg_f1_scores_weighted = [float(sum(col)) / len(col) for col in zip(*metrics['all_f1_scores']['weighted'])]

    max_val_acc = max(avg_val_accs_across_folds)
    max_val_index = avg_val_accs_across_folds.index(max_val_acc)
    loss = avg_loss_across_folds[max_val_index]
    best_epoch = max_val_index * config.eval_model_every if max_val_index != 1 else 1
    best_epoch_index = int(best_epoch / config.eval_model_every)


    best_class_f1_scores = {c: [] for c in ['long-lived', 'normal-lived', 'short-lived']}
    for fold in range(num_folds):
        xy_val_fold = val_folds_data[fold]
        y_true = [data.y.item() for data in xy_val_fold]
        y_pred = fold_predictions[fold][best_epoch_index]


        y_pred = fold_predictions[fold][best_epoch_index]
        for class_name, class_index in zip(['long-lived', 'normal-lived', 'short-lived'], [0, 1, 2]):
            class_indices = [i for i, val in enumerate(y_true) if val == class_index]
            class_y_true = [y_true[i] for i in class_indices]
            class_y_pred = [y_pred[i] for i in class_indices]
            class_f1 = f1_score(class_y_true, class_y_pred, average=None, zero_division=0)
            best_class_f1_scores[class_name].append(np.mean(class_f1))

    best_class_f1_scores_avg = {cls: np.mean(vals) for cls, vals in best_class_f1_scores.items()}
    individual_fold_performances = [row[best_epoch_index] for row in metrics['all_val_accs']]
    best_predictions_per_fold = [preds[best_epoch_index] for preds in fold_predictions]
    best_preds = [pred for preds in best_predictions_per_fold for pred in preds]
    best_preds = pd.DataFrame(best_preds)
    predicted_proportions = pd.DataFrame(best_preds.value_counts() / sum(best_preds.value_counts()), columns=['Predicted Proportions'])
    predicted_proportions.rename(index={0: "long-lived", 1: "normal-lived", 2: "short-lived"}, inplace=True)
    predicted_proportions.reset_index(inplace=True)
    predicted_proportions.columns = ['Longevity', 'Predicted Proportions']

    f = open(f"/data/bi1/results/neural/{experiment_name}/best_model_stats.txt", "a")
    f.write(f"Best Epoch {best_epoch}. Loss: {loss}. Best Average Validation Accuracy: {max_val_acc}. "
            f"Best F1-scores: Micro: {avg_f1_scores_micro[max_val_index]:.4f}, Macro: {avg_f1_scores_macro[max_val_index]:.4f}, Weighted: {avg_f1_scores_weighted[max_
val_index]:.4f}\n")

    f.write("\nClass-level F1-scores (at Best Epoch across folds):\n")
    for class_name, avg_f1 in best_class_f1_scores_avg.items():
        f.write(f"  {class_name}: {avg_f1:.4f}\n")
        print(f"  {class_name}: {avg_f1:.4f}")

    if num_folds == 5:
        df_fold_performance = pd.DataFrame(np.zeros([1, 5]), columns=[f"Fold {i}" for i in range(5)])
    elif num_folds == 10:
        df_fold_performance = pd.DataFrame(np.zeros([1, 10]), columns=[f"Fold {i}" for i in range(10)])
    df_fold_performance.iloc[0, :] = individual_fold_performances
    f.write("\nPerformance across folds:\n")
    f.write(df_fold_performance.to_string(index=False))
    f.write("\n\n")
    f.write("Overall predicted proportions:\n")
    f.write(predicted_proportions.to_string(index=False))
    f.write("\n")

def train_epoch(loader, model, opt, device):
    total_loss = 0
    model.train()
    for batch in tqdm(loader):
        opt.zero_grad()
        label = batch.y
        if (model.__class__.__name__ == "MLP"):
            x = batch.x.reshape(len(batch.y), -1).to(device)
            pred = model(x)
            label = torch.LongTensor(np.concatenate(label))
        else:
            batch.to(device)
            pred = model(batch)

        loss = model.loss(pred, label.to(device))
        loss.backward()
        opt.step()
        total_loss += loss.item() * batch.num_graphs
    total_loss /= len(loader.dataset)
    return model, opt, total_loss

def test_epoch(loader, model, device):
    model.eval()
    correct = 0
    all_labels = []
    all_predictions = []
    for data in tqdm(loader):
        with torch.no_grad():
            label = data.y
            if model.__class__.__name__ == "MLP":
                x = data.x.reshape(len(data.y), -1).to(device)
                pred = model(x)
                label = torch.LongTensor(np.concatenate(label))
            else:
                data = data.to(device)
                pred = model(data)
            label = label.to(device)
            pred = pred.argmax(dim=1)
            all_labels.extend(label.cpu().numpy())  # <EC><8B><A4><EC><A0><9C> <EB><9D><BC><EB><B2><A8> <EC>

            all_predictions.extend(pred.cpu().numpy())  # <EC><98><88><EC><B8><A1><EA><B0><92> <EC><A0><80><

            correct += pred.eq(label).sum().item()  # <EC><A0><95><ED><99><95><EB><8F><84> <EA><B3><84><EC><


    total = len(loader.dataset)
    val_acc = correct / total

    # F1-score <EA><B3><84><EC><82><B0>
    f1_scores= {
        'micro': f1_score(all_labels, all_predictions, average='micro'),
        'macro': f1_score(all_labels, all_predictions, average='macro'),
        'weighted': f1_score(all_labels, all_predictions, average='weighted')
    }

    class_report = classification_report(all_labels, all_predictions, output_dict=True)
    class_f1_scores = {
        k: v['f1-score'] for k, v in class_report.items()
        if k in ['0', '1', '2']  # <ED><81><B4><EB><9E><98><EC><8A><A4> <EB><9D><BC><EB><B2><A8><EC><9D><B4>
 
    }

    return val_acc, all_predictions, f1_scores, class_f1_scores


def load_data(config):
    dataset = DatasetFromCSV(f'{config.expression_path}',
                                 f'{config.label_path}',
                                 f'{config.age_path}',
                                 f'{config.experiments_path}',
                                 config.aging_genes_only)

    try:
        if (config.train_GNN):
import numpy as np
import torch
import os
from torch_geometric.data import DataLoader, Data
from tqdm import tqdm
from split import StratifiedGroupKFold
from pathlib import Path
from sklearn.model_selection import StratifiedShuffleSplit
import pandas as pd
from data.gene_expression_dataset import DatasetFromCSV
from data.gene_interaction_graph import StringDBGraph
from data.gene_expression_graph import GeneExpressionGraph
from sklearn.metrics import f1_score,classification_report



def cross_validate(model, X, labels, device, config, experiments, experiment_name,
                   num_folds, learning_rate, batch_size, weight_decay, epochs,
                   mixsplit=False):
    Path(f"gnn_checkpoints").mkdir(exist_ok=True)
    if mixsplit:
        stratified_k_fold = StratifiedShuffleSplit(num_folds, test_size=0.1, random_state=config.seed)
    else:
        stratified_k_fold = StratifiedGroupKFold(num_folds)

    device = device
    metrics = {
        'all_val_accs': [],
        'all_losses': [],
        'all_f1_scores': {'micro': [], 'macro': [], 'weighted': []},
        'class_metrics': {'long-lived': [], 'normal-lived': [], 'short-lived': []}
    }
    fold_predictions = []
    val_indices_list = []
    val_folds_data = [] 

    for fold, (train_index, val_index) in enumerate(stratified_k_fold.split(np.zeros(len(labels)), labels.values.reshape(-1), experiments)):
        val_indices_list.append(val_index)
        Path(f"/data/bi1/results/neural/{experiment_name}").mkdir(parents=True, exist_ok=True)
        f = open(f"/data/bi1/results/neural/{experiment_name}/fold_{fold}_training_progress.txt", "w")
        model.apply(model.weight_reset)

        if config.train_GNN:
            xy_train_fold = X[list(train_index)]
            xy_val_fold = X[list(val_index)]
        else:
            xy_train_fold = [Data(x=torch.Tensor(x.values), y=labels.loc[idx].values) for idx, x in X.iloc[list(train_index)].iterrows()]
            xy_val_fold = [Data(x=torch.Tensor(x.values), y=labels.loc[idx].values) for idx, x in X.iloc[list(val_index)].iterrows()]
        val_folds_data.append(xy_val_fold)

        opt = torch.optim.Adam(model.parameters(), lr=learning_rate, weight_decay=weight_decay)
        train_loader = DataLoader(xy_train_fold, batch_size=batch_size, shuffle=True)
        val_loader = DataLoader(xy_val_fold, batch_size=batch_size, shuffle=True)

        fold_val_accs = []
        fold_losses = []
        fold_f1_scores_micro = []
        fold_f1_scores_macro = []
        fold_f1_scores_weighted = []
        all_epoch_predictions = []
        loaded_upper_checkpoint = False

        for epoch in range(1, epochs + 1):
            upper_checkpoint = int((epoch + config.eval_model_every) / config.eval_model_every) * config.eval_model_every if epoch % config.eval_model_every != 0 and 
epoch != 1 else epoch
            upper_checkpoint_path = f"gnn_checkpoints/fold_{fold}_epoch_{upper_checkpoint}_{experiment_name}.pt.tar"
            path = f"gnn_checkpoints/fold_{fold}_epoch_{epoch}_{experiment_name}.pt.tar"

            if os.path.exists(upper_checkpoint_path):
                if not loaded_upper_checkpoint:
                    checkpoint = torch.load(upper_checkpoint_path)
                    model.load_state_dict(checkpoint['model_state_dict'])
                    opt.load_state_dict(checkpoint['optimizer_state_dict'])
                    total_loss = checkpoint['loss']
                    loaded_upper_checkpoint = True
            else:
                model, opt, total_loss = train_epoch(train_loader, model, opt, device)

            if epoch % config.eval_model_every == 0 or epoch == 1:
                if config.train_GNN and not loaded_upper_checkpoint:
                    torch.save({
                        'epoch': epoch,
                        'model_state_dict': model.state_dict(),
                        'optimizer_state_dict': opt.state_dict(),
                        'loss': total_loss,
                    }, path)
                loaded_upper_checkpoint = False

                val_acc, epoch_predictions, f1_scores, class_f1_scores = test_epoch(val_loader, model, device)
                y_true = [data.y.item() for data in xy_val_fold]
                y_pred = epoch_predictions

                for class_name, class_index in zip(['long-lived', 'normal-lived', 'short-lived'], [0, 1, 2]):
                    class_indices = [i for i, val in enumerate(y_true) if val == class_index]
                    class_y_true = [y_true[i] for i in class_indices]
                    class_y_pred = [y_pred[i] for i in class_indices]
                    class_f1 = f1_score(class_y_true, class_y_pred, average=None, zero_division=0)
                    metrics['class_metrics'][class_name].append(np.mean(class_f1))

                overall_micro_f1 = f1_score(y_true, y_pred, average='micro', zero_division=0)
                overall_macro_f1 = f1_score(y_true, y_pred, average='macro', zero_division=0)
                overall_weighted_f1 = f1_score(y_true, y_pred, average='weighted', zero_division=0)

                fold_val_accs.append(val_acc)
                fold_losses.append(total_loss)
                fold_f1_scores_micro.append(overall_micro_f1)
                fold_f1_scores_macro.append(overall_macro_f1)
                fold_f1_scores_weighted.append(overall_weighted_f1)
                all_epoch_predictions.append(epoch_predictions)

                f.write(f"Fold {fold}. Epoch {epoch}. Loss: {total_loss}. Validation accuracy: {val_acc}. \n")

        metrics['all_val_accs'].append(fold_val_accs)
        metrics['all_losses'].append(fold_losses)
        metrics['all_f1_scores']['micro'].append(fold_f1_scores_micro)
        metrics['all_f1_scores']['macro'].append(fold_f1_scores_macro)
        metrics['all_f1_scores']['weighted'].append(fold_f1_scores_weighted)
        fold_predictions.append(all_epoch_predictions)

    avg_val_accs_across_folds = [float(sum(col)) / len(col) for col in zip(*metrics['all_val_accs'])]
    avg_loss_across_folds = [float(sum(col)) / len(col) for col in zip(*metrics['all_losses'])]
    avg_f1_scores_micro = [float(sum(col)) / len(col) for col in zip(*metrics['all_f1_scores']['micro'])]
    avg_f1_scores_macro = [float(sum(col)) / len(col) for col in zip(*metrics['all_f1_scores']['macro'])]
    avg_f1_scores_weighted = [float(sum(col)) / len(col) for col in zip(*metrics['all_f1_scores']['weighted'])]

    max_val_acc = max(avg_val_accs_across_folds)
    max_val_index = avg_val_accs_across_folds.index(max_val_acc)
    loss = avg_loss_across_folds[max_val_index]
    best_epoch = max_val_index * config.eval_model_every if max_val_index != 1 else 1
    best_epoch_index = int(best_epoch / config.eval_model_every)


    best_class_f1_scores = {c: [] for c in ['long-lived', 'normal-lived', 'short-lived']}
    for fold in range(num_folds):
        xy_val_fold = val_folds_data[fold]
        y_true = [data.y.item() for data in xy_val_fold]
        y_pred = fold_predictions[fold][best_epoch_index]


        y_pred = fold_predictions[fold][best_epoch_index]
        for class_name, class_index in zip(['long-lived', 'normal-lived', 'short-lived'], [0, 1, 2]):
            class_indices = [i for i, val in enumerate(y_true) if val == class_index]
            class_y_true = [y_true[i] for i in class_indices]
            class_y_pred = [y_pred[i] for i in class_indices]
            class_f1 = f1_score(class_y_true, class_y_pred, average=None, zero_division=0)
            best_class_f1_scores[class_name].append(np.mean(class_f1))

    best_class_f1_scores_avg = {cls: np.mean(vals) for cls, vals in best_class_f1_scores.items()}
    individual_fold_performances = [row[best_epoch_index] for row in metrics['all_val_accs']]
    best_predictions_per_fold = [preds[best_epoch_index] for preds in fold_predictions]
    best_preds = [pred for preds in best_predictions_per_fold for pred in preds]
    best_preds = pd.DataFrame(best_preds)
    predicted_proportions = pd.DataFrame(best_preds.value_counts() / sum(best_preds.value_counts()), columns=['Predicted Proportions'])
    predicted_proportions.rename(index={0: "long-lived", 1: "normal-lived", 2: "short-lived"}, inplace=True)
    predicted_proportions.reset_index(inplace=True)
    predicted_proportions.columns = ['Longevity', 'Predicted Proportions']

    f = open(f"/data/bi1/results/neural/{experiment_name}/best_model_stats.txt", "a")
    f.write(f"Best Epoch {best_epoch}. Loss: {loss}. Best Average Validation Accuracy: {max_val_acc}. "
            f"Best F1-scores: Micro: {avg_f1_scores_micro[max_val_index]:.4f}, Macro: {avg_f1_scores_macro[max_val_index]:.4f}, Weighted: {avg_f1_scores_weighted[max_
val_index]:.4f}\n")

    f.write("\nClass-level F1-scores (at Best Epoch across folds):\n")
    for class_name, avg_f1 in best_class_f1_scores_avg.items():
        f.write(f"  {class_name}: {avg_f1:.4f}\n")
        print(f"  {class_name}: {avg_f1:.4f}")

    if num_folds == 5:
        df_fold_performance = pd.DataFrame(np.zeros([1, 5]), columns=[f"Fold {i}" for i in range(5)])
    elif num_folds == 10:
        df_fold_performance = pd.DataFrame(np.zeros([1, 10]), columns=[f"Fold {i}" for i in range(10)])
    df_fold_performance.iloc[0, :] = individual_fold_performances
    f.write("\nPerformance across folds:\n")
    f.write(df_fold_performance.to_string(index=False))
    f.write("\n\n")
    f.write("Overall predicted proportions:\n")
    f.write(predicted_proportions.to_string(index=False))
    f.write("\n")

def train_epoch(loader, model, opt, device):
    total_loss = 0
    model.train()
    for batch in tqdm(loader):
        opt.zero_grad()
        label = batch.y
        if (model.__class__.__name__ == "MLP"):
            x = batch.x.reshape(len(batch.y), -1).to(device)
            pred = model(x)
            label = torch.LongTensor(np.concatenate(label))
        else:
            batch.to(device)
            pred = model(batch)

        loss = model.loss(pred, label.to(device))
        loss.backward()
        opt.step()
        total_loss += loss.item() * batch.num_graphs
    total_loss /= len(loader.dataset)
    return model, opt, total_loss

def test_epoch(loader, model, device):
    model.eval()
    correct = 0
    all_labels = []
    all_predictions = []
    for data in tqdm(loader):
        with torch.no_grad():
            label = data.y
            if model.__class__.__name__ == "MLP":
                x = data.x.reshape(len(data.y), -1).to(device)
                pred = model(x)
                label = torch.LongTensor(np.concatenate(label))
            else:
                data = data.to(device)
                pred = model(data)
            label = label.to(device)
            pred = pred.argmax(dim=1)
            all_labels.extend(label.cpu().numpy())  # <EC><8B><A4><EC><A0><9C> <EB><9D><BC><EB><B2><A8> <EC>

            all_predictions.extend(pred.cpu().numpy())  # <EC><98><88><EC><B8><A1><EA><B0><92> <EC><A0><80><

            correct += pred.eq(label).sum().item()  # <EC><A0><95><ED><99><95><EB><8F><84> <EA><B3><84><EC><


    total = len(loader.dataset)
    val_acc = correct / total

    # F1-score <EA><B3><84><EC><82><B0>
    f1_scores= {
        'micro': f1_score(all_labels, all_predictions, average='micro'),
        'macro': f1_score(all_labels, all_predictions, average='macro'),
        'weighted': f1_score(all_labels, all_predictions, average='weighted')
    }

    class_report = classification_report(all_labels, all_predictions, output_dict=True)
    class_f1_scores = {
        k: v['f1-score'] for k, v in class_report.items()
        if k in ['0', '1', '2']  # <ED><81><B4><EB><9E><98><EC><8A><A4> <EB><9D><BC><EB><B2><A8><EC><9D><B4>
 
    }

    return val_acc, all_predictions, f1_scores, class_f1_scores


def load_data(config):
    dataset = DatasetFromCSV(f'{config.expression_path}',
                                 f'{config.label_path}',
                                 f'{config.age_path}',
                                 f'{config.experiments_path}',
                                 config.aging_genes_only)

    try:
        if (config.train_GNN):
            gene_interaction_graph = StringDBGraph()
            gene_expression_graphs = GeneExpressionGraph(gene_interaction_graph, dataset.df, config, dataset
.labels,
                                                         dataset.age)
            pre_expression_merge_graph = gene_expression_graphs.gene_exp_graph_pre_exp_values_merge
            return dataset, gene_expression_graphs, pre_expression_merge_graph
        else:
            return dataset
    except:
        return dataset

def load_train_test_datasets(dataset, gene_expression_graphs=None):
    cv = StratifiedGroupKFold(10)
    labels = dataset.labels.to_numpy().reshape(-1)

    for idx, (train_idxs, test_idxs) in enumerate(cv.split(dataset.df.to_numpy(), labels,
                                                           dataset.experiments)):
        if (idx == 1):
            break
        if (gene_expression_graphs):
            X_train = gene_expression_graphs[train_idxs]
            X_test = gene_expression_graphs[test_idxs]
        else:
            age = dataset.age
            X = dataset.df.merge(age, left_index=True, right_index=True)
            X_train = X.iloc[train_idxs]
            X_test = X.iloc[test_idxs]

        labels_train = dataset.labels.iloc[train_idxs]
        experiments_train = dataset.experiments.iloc[train_idxs]
        labels_test = dataset.labels.iloc[test_idxs]
        experiments_test = dataset.experiments.iloc[test_idxs]

    return [(X_train, labels_train, experiments_train),
            (X_test, labels_test, experiments_test)]

def get_num_genes(dataset):
    return dataset.df.shape[1]
