from pathlib import Path
import pandas as pd
import numpy as np
import scanpy as sc

# Locate the project root based on this file's location
repo_root = Path(__file__).resolve().parents[1]

# Key folders (as Path objects)
data_dir = repo_root / "data"
src_dir = data_dir / "GSE138852"   # where CSVs are
proc_dir = data_dir / "processed"   # where .h5ad will be saved
fig_dir = repo_root / "figures"
res_dir = repo_root / "results"

# Create folders to write to
for d in (proc_dir, fig_dir, res_dir):
    d.mkdir(parents=True, exist_ok=True)

print("Repo root:", repo_root)
print("Output folders ready:", proc_dir, fig_dir, res_dir)
