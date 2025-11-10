# Alzheimer's scRNA-seq Analysis

This repository contains a bioinformatics project analyzing single-cell RNA sequencing (scRNA-seq) data from the entorhinal cortex of Alzheimer's disease (AD) patients. The goal is to process and explore transcriptional changes across individual cells and patient samples, using standard single-cell workflows in Python.

## Objectives

- Load and explore scRNA-seq gene expression and metadata
- Perform quality control, normalization, and dimensionality reduction
- Identify cell types and marker genes
- Compare expression profiles across patient groups
- Visualize results using UMAP, heatmaps, and other techniques

## Data

- Dataset: GSE138852 (processed counts and covariates)
- Platform: Illumina NextSeq 500
- Tissue: Human entorhinal cortex from aged individuals with AD

Raw data was excluded from version control due to file size limits. Metadata and scripts are included for reproducibility.

## Tools and Packages

- **Programming languages**: Python, Bash
- **Python libraries for data science**: pandas, numpy, matplotlib, seaborn
- **Single-cell and genomics analysis**: scanpy, anndata, celltypist, scvi-tools
- **Statistical modeling**: statsmodels, lifelines
- **Bioinformatics tools**: samtools, bcftools, bedtools, bwa, BLAST (used outside Python)
- **File formats handled**: FASTQ, BAM, VCF, BED, BigWig, H5AD, CSV/TSV
- **Databases**: PostgreSQL, SQL for metadata storage and queries
- **Version control**: Git, GitHub
- **Workflow management**: Snakemake (planned)
- **Visualization**: matplotlib, seaborn, scanpy plotting utilities

## Skills Demonstrated

This project showcases:

- **Programming:** Python (pandas, numpy, scanpy, anndata, matplotlib, seaborn) and Bash scripting for data preprocessing
- **Data Science:** Data cleaning, wrangling, and exploratory data analysis on high-dimensional biological data
- **Bioinformatics:** Single-cell RNA-seq preprocessing, QC, normalization, clustering, and marker gene identification
- **NGS Tools & Formats:** Familiarity with common genomics file formats (FASTQ, BAM, VCF, H5AD) and tools such as samtools bcftools, bedtools
- **Statistics & Machine Learning:** Dimensionality reduction (PCA, UMAP), clustering (Leiden), differential expression analysis
- **Reproducibility & Workflow Management:** Version control with Git/GitHub, clear directory structure, reproducible Python scripts and Jupyter notebooks

## Project Structure

```

scrnaseq/
├── config/        # Configuration files for workflows and analysis
├── data/          # Project data directory
│   ├── raw/       # Raw input datasets (e.g., FASTQ, REDCap exports)
│   ├── ref/       # Reference data (e.g., genome fasta, annotations)
│   ├── processed/ # Processed datasets (e.g., AnnData, intermediate CSVs)
│   └── GSE138852/ # Alzheimer’s counts and covariates (from GEO)
├── figures/       # Generated plots (excluded from Git)
├── notebooks/     # Jupyter notebooks for interactive analysis
├── reports/       # Generated HTML/PDF reports (excluded from Git)
├── results/       # Analysis outputs and intermediate files (excluded from Git)
├── rules/         # Workflow engine rules (e.g., Snakemake)
├── scripts/       # Python scripts for loading, processing, and visualization
├── sql/           # Database schema and query files
└── README.md      # Project overview
```

## Status

Initial data loading in Python complete. QC and analysis notebooks in development.