#!/usr/bin/env bash
# clean.sh — remove generated and temporary files from the project

set -e  # exit if any command fails

echo "Cleaning processed data, figures, results, and temp files..."

rm -rf data/processed/* figures/* results/* reports/* .snakemake/ 2>/dev/null || true
find . -name ".DS_Store" -delete
find . -name ".ipynb_checkpoints" -type d -exec rm -rf {} +

echo "Cleanup complete."

