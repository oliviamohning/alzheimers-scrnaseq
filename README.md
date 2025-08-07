# Alzheimer's scRNA-seq Analysis

This repository contains a bioinformatics project analyzing single-cell RNA sequencing (scRNA-seq) data from the entorhinal cortex of Alzheimer's disease (AD) patients. The goal is to process and explore transcriptional changes across individual cells and patient samples, using standard single-cell workflows in R (Seurat).

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

- R / RStudio
- Seurat
- tidyverse (dplyr, ggplot2, etc.)
- Future additions: DESeq2, clusterProfiler, etc.

## Project Structure

```
scrnaseq/
├── data/       # Processed datasets (counts, covariates, saved Seurat objects)
├── scripts/    # R scripts for loading, processing, and visualizing data
├── results/    # Outputs and intermediate files (excluded from Git)
├── figures/    # Generated plots (excluded from Git)
└── README.md   # Project overview
```

## Status

Initial data loading and object creation complete. Analysis scripts in development.
