# Load required libraries
library(Seurat)
library(tidyverse)

# Define file paths
counts_path <- "/Users/oliviamohning/Documents/ds-portfolio/scrnaseq/data/GSE138852/GSE138852_counts.csv"
meta_path   <- "/Users/oliviamohning/Documents/ds-portfolio/scrnaseq/data/GSE138852/GSE138852_covariates.csv"

# Load count matrix
counts <- read.csv(counts_path, row.names = 1)

# Load metadata
metadata <- read.csv(meta_path, row.names = 1)

# Check dimensions
dim(counts)     # genes x cells
dim(metadata)   # cells x metadata columns

# Create Seurat object
alz_obj <- CreateSeuratObject(counts = counts, meta.data = metadata)

# Save the Seurat object to file
saveRDS(alz_obj, file = "/Users/oliviamohning/Documents/ds-portfolio/scrnaseq/data/GSE138852/alz_seurat.rds")
