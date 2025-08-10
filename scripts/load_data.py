# Silence some warnings
import warnings
warnings.filterwarnings("ignore", category=UserWarning, module=r"louvain")

# Importing libraries
from pathlib import Path
import pandas as pd
import numpy as np
import scanpy as sc

# Locate the project root based on this file's location
repo_root = Path(__file__).resolve().parents[1]

# Create key folders
data_dir = repo_root / "data"
src_dir = data_dir / "GSE138852"   # where CSVs are
proc_dir = data_dir / "processed"   # where .h5ad will be saved
fig_dir = repo_root / "figures"
res_dir = repo_root / "results"

# Create folders to write to
for d in (proc_dir, fig_dir, res_dir):
    d.mkdir(parents=True, exist_ok=True)

# Input files
counts_path = src_dir / "GSE138852_counts.csv"
meta_path   = src_dir / "GSE138852_covariates.csv"

# Existence checks
if not counts_path.exists():
    raise FileNotFoundError(f"Missing file: {counts_path}")
if not meta_path.exists():
    raise FileNotFoundError(f"Missing file: {meta_path}")

# Read inputs
counts = pd.read_csv(counts_path, index_col=0)
meta   = pd.read_csv(meta_path,   index_col=0)

# Transpose counts so rows are cells and columns are genes
expr_matrix = counts.T

# Align metadata with expression matrix
common_cells = expr_matrix.index.intersection(meta.index)
expr_matrix = expr_matrix.loc[common_cells].copy()
meta = meta.loc[common_cells].copy()

# Create AnnData object
adata = sc.AnnData(X=expr_matrix)
adata.obs = meta
adata.var.index = expr_matrix.columns
adata.var.index.name = "gene"
adata.obs.index.name = "cell"

# Save AnnData object
h5ad_path = proc_dir / "alz_scanpy_raw.h5ad"
adata.write(h5ad_path, compression="gzip")
print(f"Saved AnnData to: {h5ad_path}")
