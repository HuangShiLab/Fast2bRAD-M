#!/bin/bash
#SBATCH --job-name=ecc_downstream
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=04:00:00
#SBATCH --output=ecc_downstream_%j.log
set -euo pipefail

# Activate the holo2bRAD Python venv
source /lustre1/g/aos_shihuang/holo2bRAD/.venv/bin/activate

SCRIPT_DIR="/lustre1/g/aos_shihuang/holo2bRAD/caries_pilot"
METADATA="${SCRIPT_DIR}/metadata.tsv"
HOLO_DIR="/lustre1/g/aos_shihuang/holo2bRAD/caries_pilot/holo/holo"
OUTDIR="/lustre1/g/aos_shihuang/holo2bRAD/caries_pilot/results"

mkdir -p "${OUTDIR}"

echo "Building ECC feature tables..."
python3 "${SCRIPT_DIR}/build_ecc_features.py" \
    --metadata "${METADATA}" \
    --holo-dir "${HOLO_DIR}" \
    --enzyme BcgI \
    --output "${OUTDIR}/feature_tables"

echo "Running classification benchmark..."
python3 "${SCRIPT_DIR}/classify_ecc.py" \
    --feature-dir "${OUTDIR}/feature_tables" \
    --output "${OUTDIR}/classification"

echo "Done. Results in ${OUTDIR}/"
