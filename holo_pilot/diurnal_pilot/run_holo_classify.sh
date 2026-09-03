#!/bin/bash
#SBATCH --job-name=holo_diurnal
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=24:00:00
#SBATCH --output=holo_diurnal_%j.log
set -euo pipefail

# Usage: run_holo_classify.sh <metadata.tsv> <enzyme> <output_dir> <threads>
METADATA="${1:?metadata.tsv required}"
ENZYME="${2:?enzyme required, e.g. BcgI}"
OUTDIR="${3:?output directory required}"
THREADS="${4:-16}"

# ----- HPC PATHS -----
F2BRAD_HOLO="/lustre1/g/aos_shihuang/Fast2bRAD-M/target/release/f2brad-holo"
HOST_DB="/lustre1/g/aos_shihuang/holo2bRAD/host_db/${ENZYME}.host_db.tsv"
MICROBE_DB="/lustre1/g/aos_shihuang/Fast2bRAD-M/db/02_db_qual/${ENZYME}.species.iibdb"
MICROBE_MASK="/lustre1/g/aos_shihuang/holo2bRAD/cross_results/${ENZYME}/microbe_mask.${ENZYME}.2.txt"
MICROBE_DB_DIR="/lustre1/g/aos_shihuang/Fast2bRAD-M/db/02_db_qual"

mkdir -p "${OUTDIR}"

SAMPLE_LIST="${OUTDIR}/sample_list.tsv"
awk -F'\t' 'NR>1 && $4!="" {print $1"\t"$4}' "${METADATA}" > "${SAMPLE_LIST}"

echo "Running f2brad-holo classify for $(wc -l < "${SAMPLE_LIST}") samples"
"${F2BRAD_HOLO}" classify \
    -d "${HOST_DB}" \
    -m "${MICROBE_DB}" \
    --microbe-mask "${MICROBE_MASK}" \
    --microbe-db-dir "${MICROBE_DB_DIR}" \
    --exclude-human \
    -l "${SAMPLE_LIST}" \
    -s "${ENZYME}" \
    -o "${OUTDIR}/holo" \
    -j "${THREADS}"

echo "Done. Outputs in ${OUTDIR}/holo/<sample_id>/"
