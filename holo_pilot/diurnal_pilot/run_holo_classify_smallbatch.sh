#!/bin/bash
#SBATCH --job-name=holo_diurnal_small
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=04:00:00
#SBATCH --output=holo_diurnal_small_%j.log
set -euo pipefail

# Small-batch pilot: first 3 diurnal samples
METADATA="/lustre1/g/aos_shihuang/holo2bRAD/diurnal_pilot/metadata.tsv"
ENZYME="BcgI"
OUTDIR="/lustre1/g/aos_shihuang/holo2bRAD/diurnal_pilot/smallbatch"
THREADS="${1:-8}"

# ----- HPC PATHS -----
F2BRAD_HOLO="/lustre1/g/aos_shihuang/Fast2bRAD-M/target/release/f2brad-holo"
HOST_DB="/lustre1/g/aos_shihuang/holo2bRAD/host_db/${ENZYME}.host_db.tsv"
MICROBE_DB="/lustre1/g/aos_shihuang/Fast2bRAD-M/db/02_db_qual/${ENZYME}.species.iibdb"
MICROBE_MASK="/lustre1/g/aos_shihuang/holo2bRAD/cross_results/${ENZYME}/microbe_mask.${ENZYME}.2.txt"
MICROBE_DB_DIR="/lustre1/g/aos_shihuang/Fast2bRAD-M/db/02_db_qual"

mkdir -p "${OUTDIR}"

# First 3 samples, single-end (2 columns)
SAMPLE_LIST="${OUTDIR}/sample_list.tsv"
awk -F'\t' 'NR>1 && NR<=4 && $4!="" {print $1"\t"$4}' "${METADATA}" > "${SAMPLE_LIST}"

echo "Running small-batch f2brad-holo classify for $(wc -l < "${SAMPLE_LIST}") samples"
"${F2BRAD_HOLO}" classify \
    -d "${HOST_DB}" \
    -m "${MICROBE_DB}" \
    --microbe-mask "${MICROBE_MASK}" \
    --microbe-db-dir "${MICROBE_DB_DIR}" \
    -l "${SAMPLE_LIST}" \
    -s "${ENZYME}" \
    -o "${OUTDIR}/holo" \
    -j "${THREADS}"

echo "Done. Outputs in ${OUTDIR}/holo/<sample_id>/"
