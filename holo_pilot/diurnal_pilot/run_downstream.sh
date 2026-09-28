#!/bin/bash
#SBATCH --job-name=diurnal_downstream
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=04:00:00
#SBATCH --output=diurnal_downstream_%j.log
set -euo pipefail

# Activate the holo2bRAD Python venv and R conda env
source /lustre1/g/aos_shihuang/holo2bRAD/.venv/bin/activate
source /group/aos_shihuang/conda/etc/profile.d/conda.sh
conda activate R

SCRIPT_DIR="/lustre1/g/aos_shihuang/holo2bRAD/diurnal_pilot"
METADATA="${SCRIPT_DIR}/metadata.tsv"
HOLO_DIR="/lustre1/g/aos_shihuang/holo2bRAD/diurnal_pilot/holo/holo"
OUTDIR="/lustre1/g/aos_shihuang/holo2bRAD/diurnal_pilot/results"

mkdir -p "${OUTDIR}"

echo "Building diurnal feature tables..."
python3 "${SCRIPT_DIR}/build_diurnal_features.py" \
    --metadata "${METADATA}" \
    --holo-dir "${HOLO_DIR}" \
    --output "${OUTDIR}/feature_tables"

echo "Running diurnal host-microbe analysis..."
Rscript "${SCRIPT_DIR}/diurnal_analysis.R" \
    "${OUTDIR}/feature_tables" \
    "${OUTDIR}/diurnal_analysis"

echo "Done. Results in ${OUTDIR}/"
