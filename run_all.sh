#!/usr/bin/env bash
# Runs the whole pipeline in order and stops at the first failure.
set -euo pipefail

for s in 01_prepare 02_features 03_cluster 04_model 05_figures; do
  echo "=== $s ==="
  Rscript "analysis/$s.R"
done

echo "=== python cross-checks ==="
python3 python/cluster_check.py
python3 python/model_check.py

echo
echo "Done. Tables in outputs/, figures in figures/."
