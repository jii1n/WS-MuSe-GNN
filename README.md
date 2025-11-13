# model(제목 수정)

It contains the metadata, scripts and  model-ready demo dataset 
used in the study : < 논문제목 > 


The full processed gene expression matrix
('getmm_combat_seq_no_outliers_and_singles_gene_expression.csv')
is not included due to file size restrictions.  
However, metadata, scripts, and execution-ready demo dataset are 
provided, allowing users   
to:
- Run the machine learning pipeline immediately.
- Reproduce the full dataset from public SRA resources if desired.


## 📁 Data

```
.
├── accession_list.txt               # SRA Run IDs to download
├── age.csv                          # Sample age metadata
├── labels.csv                       # Lifespan class labels (long/normal/short)
├── sra_to_bioproject.csv            # Run → BioProject mapping
├── all_bioproject.txt               # All BioProject IDs retrieved
├── genage_all_genes.txt             # Aging gene list
├── keywords_metadata.p              # GEO keyword-based metadata
├── genage_genes_metadata.p          # GenAge-based metadata
├── manually_fetched_metadata.p      # Additional manually curated metadata
├── get_accession_list.py            # Script for metadata aggregation
├── combine_data.R                   # Combine all Kallisto outputs
├── getmm_and_combat_seq.R           # GeTMM + ComBat-seq normalization pipeline
├── create_dummy.ipynb               # Script that generated the demo dataset
├── demo_data/
│   ├── gene_exp.csv
│   ├── labels.csv
│   ├── age.csv
│   └── sra_to_bioproject.csv
```


## 🧪 method 1. Execution-ready Demo Dataset 

To run the ML models without downloading FASTQ files, 

Example
```bash
python3 cv_neural.py \
  --expression_path gene_exp.csv \
  --label_path labels.csv \
  --age_path age.csv \
  --experiments_path sra_to_bioproject.csv
```
This dataset preserves the exact data structure expected by the ML pipeline.


## 🔧 method 2. How to Reproduce the Full Expression Matrix

The full expression matrix can be regenerated from public SRA data using the  
steps below. All required metadata and processing scripts are included.

### Step 1 — Download FASTQ files
Use IDs from:
- `accession_list.txt`
- `sra_to_bioproject.csv`
- metadata files (*.p)
  
### Step 2 — Quantify reads using Kallisto

### Step 3 — Merge all quantification outputs
```
Rscript combine_data.R
```

Produces:
```
combined_studies_raw_counts.Rdata
```
  
### Step 4 — Normalize + batch correct

```
Rscript getmm_and_combat_seq.R  
```

Produces:
```
getmm_combat_seq_no_outliers_and_singles_gene_expression.csv
combat_seq_gene_expression_no_outliers_and_singles.Rdata
etc
```

## 🚀 3. Running the Machine Learning Model

A ready-to-use shell script (`mlp_runs.sh`, `focla_loss_mlp.sh`, `gnn_runs.sh`, `ws_muse_gnn_runs.sh`) is provided to execute the ML
pipeline.

To run:

```
bash mlp_runs.sh
```

Example (inside the script):

```bash
python3 cv_neural.py \
    --expression_path demo_data/gene_exp_dummy.csv \
    --label_path demo_data/labels_dummy.csv \
    --age_path demo_data/age_dummy.csv \
    --experiments_path demo_data/sra_to_bioproject_dummy.csv \
    --mlp_hidden_dim 1024 \
    --learning_rate 0.0001 \
    --num_mlp_layers 3 \
    --train_MLP \
    --dropout 0.5
```
