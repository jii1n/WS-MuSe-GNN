# model(제목 수정)

It contains the metadata, scripts and  model-ready demo dataset used in < 논문제목 > study


The full processed gene expression matrix
(`getmm_combat_seq_no_outliers_and_singles_gene_expression.csv`)
is not included due to file size restrictions.  
However, **metadata, scripts, and execution-ready demo dataset** are provided, allowing users   
to:
- Run the machine learning pipeline immediately.
- Reproduce the full dataset from public SRA resources if desired.


## 📁 data

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


## 🧪 1. Execution-ready Demo Dataset 

To run the ML models without downloading large FASTQ files, 

```bash
ex)
python3 cv_neural.py \
  --expression_path gene_exp.csv \
  --label_path labels.csv \
  --age_path age.csv \
  --experiments_path sra_to_bioproject.csv

This dataset preserves the exact data structure expected by the ML pipeline.




