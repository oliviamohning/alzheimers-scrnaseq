# Alzheimer's Disease snRNA-seq Analysis (Entorhinal Cortex)

This repository contains a bioinformatics project analyzing **single-nucleus RNA sequencing (snRNA-seq)** data from the entorhinal cortex of individuals with Alzheimer’s disease (AD) and control samples. The goal is to process and explore transcriptional changes across individual nuclei and patient groups using Python-based data science workflows.

## Objectives
- Load and explore snRNA-seq gene expression counts and metadata  
- Perform basic data wrangling and quality inspection  
- Compare cell-type composition between Alzheimer’s and control samples  
- Visualize gene expression trends for key genes (APOE, CLU, GFAP)  
- Generate reproducible plots and summary reports

## Data
- **Dataset:** GSE138852 (processed counts and covariates)  
- **Platform:** Illumina NextSeq 500  
- **Tissue:** Human entorhinal cortex from aged individuals with and without Alzheimer’s disease  

Raw data was excluded from version control due to file size limits. Metadata and processed CSVs are included for reproducibility.

## Tools and Packages
- **Programming languages:** Python, Bash  
- **Python libraries:** pandas, numpy, seaborn, matplotlib  
- **Version control:** Git, GitHub  
- **Visualization:** seaborn, matplotlib  

## Project Structure
scrnaseq/
├── config/ # Configuration files for workflows and analysis
├── data/ # Project data directory
│ ├── raw/ # Raw input datasets (e.g., FASTQ, REDCap exports)
│ ├── ref/ # Reference data (e.g., genome fasta, annotations)
│ ├── processed/ # Processed datasets (e.g., CSVs, AnnData files)
│ └── GSE138852/ # Alzheimer’s counts and covariates (from GEO)
├── figures/ # Generated plots (included in Git)
├── notebooks/ # Jupyter notebooks for interactive analysis
├── reports/ # Markdown or HTML reports summarizing findings
├── results/ # Analysis outputs and intermediate files
├── scripts/ # Python scripts for loading, processing, and visualization
└── README.md # Project overview


## Status
Initial data ingestion, cleaning, and exploratory analysis completed.  
Visualizations and summary report available in `/reports/summary.md`.
