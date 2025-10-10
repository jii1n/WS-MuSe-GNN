import torch
import torch.nn as nn
import torch.nn.functional as F
from torch_geometric.nn import TransformerConv, global_sort_pool
from torch_sparse import SparseTensor
from torch.nn import BatchNorm1d
from neural_models.mlp import MLP


class MuSeGNN(nn.Module):
    def __init__(self, input_dim, backbone_channels, output_dim, num_backbone_layers,
                 concat_input_graph, num_nodes, edge_attr, edge_index, k,
                 mlp_hidden_dim, mlp_num_layers, dropout, device='cpu', heads=4):
        super(MuSeGNN, self).__init__()
        self.device = torch.device(device)
        self.input_dim = input_dim
        self.num_nodes = num_nodes
        self.k = k
        self.backbone_channels = backbone_channels
        self.num_backbone_layers = num_backbone_layers
        self.dropout = dropout
        self.concat_input_graph = concat_input_graph
        self.heads = heads
        self.output_dim = output_dim

        # Convert edge_attr and edge_index to SparseTensor
        edge_attr = torch.as_tensor(edge_attr, dtype=torch.float32, device=self.device)
        self.edge_index_sp = SparseTensor(
            row=edge_index[0].to(self.device),
            col=edge_index[1].to(self.device),
            value=edge_attr,
            sparse_sizes=(num_nodes, num_nodes)
        ).to(self.device)

        # Weight sharing
        self.shared_conv = TransformerConv(
            input_dim, backbone_channels, heads=self.heads, concat=False, dropout=dropout
        )
        
        self.convs = nn.ModuleList([self.shared_conv for _ in range(num_backbone_layers)])
        self.bns = nn.ModuleList([BatchNorm1d(backbone_channels) for _ in range(num_backbone_layers)])

        # Calculate MLP input dimension
        if self.concat_input_graph:
            mlp_input_dim =  k * (num_backbone_layers * backbone_channels * heads + input_dim) + 1
        else:
            mlp_input_dim = k * num_backbone_layers * backbone_channels * heads + 1

       #print(f"Calculated MLP input dimension: {mlp_input_dim}")

       # Initialize MLP
        self.post_mp = MLP(mlp_input_dim, mlp_hidden_dim, output_dim, mlp_num_layers, dropout)

    def build_conv_model(self, input_dim, hidden_dim):
        """Build TransformerConv for MuSeGNN."""
        return TransformerConv(input_dim, hidden_dim, heads=self.heads, concat=False, dropout=self.dropout)

    def forward(self, data):
        x, batch = data.x, data.batch
        edge_index = torch.stack(self.edge_index_sp.coo()[:2]).long()
        batch_size = data.y.shape[0]

        # Reshape x to [batch_size * num_nodes, input_dim]
        x = x.reshape(batch_size * self.num_nodes, self.input_dim)
        feats_at_different_layers = []

        h = x
        for i in range(self.num_backbone_layers):
            h = self.convs[i](h, edge_index)     
            h = self.bns[i](h)                   
            h = F.relu(h)
            h = F.dropout(h, p=self.dropout, training=self.training)
            feats_at_different_layers.append(h)   

         # Concatenate features across layers
        if self.concat_input_graph:
            x = torch.cat(feats_at_different_layers, dim=-1)   
        else:
            x = feats_at_different_layers[-1]                

        # Global pooling 
        if self.k != self.num_nodes:
            x = global_sort_pool(x, batch, self.k)             
        else:
            x = x.reshape(batch_size, self.num_nodes, -1).mean(dim=1) 

        x = x.reshape(batch_size, -1)
        age = data.age.reshape(batch_size, 1).to(self.device)
        x = torch.cat((x, age), dim=1)                    


        expected_dim = self.post_mp[0].in_features
        if x.size(1) != expected_dim:
            # print(f"Adjusting MLP input dim from {x.size(1)} to {expected_dim}")
            delta = expected_dim - x.size(1)
            if delta > 0:  # padding
                x = torch.cat([x, torch.zeros(x.size(0), delta, device=x.device)], dim=1)
            elif delta < 0:  # trimming
                x = x[:, :expected_dim]

        x = self.post_mp(x)                   
        x = F.log_softmax(x, dim=-1)         
        return x

    def model_reset(self):
        for layer in self.children():
            if hasattr(layer, 'reset_parameters'):
                layer.reset_parameters()

    def weight_reset(self, m):
        if isinstance(m, nn.Linear):
            m.reset_parameters()

    def loss(self, pred, label):
        return F.nll_loss(pred, label)
