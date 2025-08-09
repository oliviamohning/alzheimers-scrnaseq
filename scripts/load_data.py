# scripts/load_data.py
# Purpose: load counts + metadata from CSV, build an AnnData object for Scanpy, and save it.

import os
import pandas as pd
import scanpy as sc

# 1) Define project paths (adjust only if you moved the folder)
proj_dir = "/Users/oliviamohning/Documents/ds-portfolio/scrnaseq"
data_dir = os.path.join(proj_dir, "data", "GSE138852")
fig_dir  = os.path.join(proj_dir, "figures")
res_dir  = os.path.join(proj_dir, "results")

os.makedirs(fig_dir, exist_ok=True)   # ensure figures/ exists
os.makedirs(res_dir, exist_ok=True)   # ensure results/ exists

# 2) File paths to CSVs you already downloaded
counts_path = os.path.join(data_dir, "GSE138852_counts.csv")
meta_path   = os.path.join(data_dir, "GSE138852_covariates.csv")

# 3) Read CSVs into pandas DataFrames
#    counts: rows=genes, cols=cells (from your R step); we will transpose for Scanpy
counts_df = pd.read_csv(counts_path, index_col=0)
meta_df   = pd.read_csv(meta_path,   index_col=0)

print("Loaded counts:", counts_df.shape, " (genes x cells)")
print("Loaded metadata:", meta_df.shape, " (cells x columns)")

# 4) Transpose counts for Scanpy (Scanpy expects cells x genes)
X = counts_df.T
print("Transposed counts:", X.shape, " (cells x genes)")

# 5) Align metadata to the cells in X
#    Ensure the metadata index matches X.index (cell IDs)
missing_meta = X.index.difference(meta_df.index)
missing_cells = meta_df.index.difference(X.index)
print("Cells in counts not in metadata:", len(missing_meta))
print("Cells in metadata not in counts:", len(missing_cells))

#    If there are mismatches, restrict both to the intersection
common_cells = X.index.intersection(meta_df.index)
X = X.loc[common_cells].copy()
meta_aligned = meta_df.loc[common_cells].copy()

print("After aligning:", X.shape, meta_aligned.shape)

# 6) Build AnnData (cells in .obs, genes in .var)
adata = sc.AnnData(X=X)
adata.obs = meta_aligned
adata.var.index.name = "gene"
adata.obs.index.name = "cell"

# 7) Basic sanity checks and lightweight summary
print(adata)
print("Observed columns:", list(adata.obs.columns)[:10])

# 8) Save AnnData to disk (large file; keep out of Git)
h5ad_path = os.path.join(data_dir, "alz_scanpy_raw.h5ad")
adata.write(h5ad_path, compression="gzip")
print("Saved AnnData to:", h5ad_path)
