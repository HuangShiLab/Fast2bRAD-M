#!/bin/bash
set -e

HPC_USER=shihuang
HPC_HOST=hpc2021.hku.hk
HPC_DIR=/home/shihuang/fast2brad_holo/cross_results_rep
LOCAL_DIR="$(cd "$(dirname "$0")/.." && pwd)/cross_results"

mkdir -p "$LOCAL_DIR"

echo "Rsyncing cross results from HPC..."
rsync -avz --progress "${HPC_USER}@${HPC_HOST}:${HPC_DIR}/" "$LOCAL_DIR/"

echo "Running summary..."
python3 "$(dirname "$0")/summarize_cross.py"

echo "Done. Results in $LOCAL_DIR and results/cross_collision_summary.md"
