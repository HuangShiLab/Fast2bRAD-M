#!/bin/bash
#SBATCH --job-name=scallop_classify
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=04:00:00
#SBATCH --output=/lustre1/g/aos_shihuang/holo2bRAD/scallop_poc/results/classify_%j.log
set -euo pipefail

REPO=/lustre1/g/aos_shihuang/Fast2bRAD-M
WORKDIR=/lustre1/g/aos_shihuang/holo2bRAD/scallop_poc
HOST_DB=$WORKDIR/db/scallop_BsaXI_host_db.tsv
MICROBE_DB=$WORKDIR/db/BsaXI_example_microbe_db/BsaXI.species.quant.iibdb
MICROBE_DIR=$WORKDIR/db/BsaXI_example_microbe_db
SAMPLES=$WORKDIR/reads/samples.tsv
OUT=$WORKDIR/results/holo

mkdir -p "$OUT"

$REPO/target/release/f2brad-holo classify \
    --host-db "$HOST_DB" \
    --microbe-db "$MICROBE_DB" \
    --microbe-db-dir "$MICROBE_DIR" \
    --site BsaXI \
    --sample-list "$SAMPLES" \
    --output "$OUT" \
    --taxonomy species \
    --host-max-mismatch 2 \
    --min-depth 4 \
    -j 16

echo "Done. Results in $OUT"
