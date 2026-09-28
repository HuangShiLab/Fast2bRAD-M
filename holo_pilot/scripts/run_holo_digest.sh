#!/bin/bash
#SBATCH --job-name=holo_digest
#SBATCH --partition=amd
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=04:00:00
#SBATCH --output=/home/shihuang/fast2brad_holo/holo_digest_%j.out
#SBATCH --error=/home/shihuang/fast2brad_holo/holo_digest_%j.err

set -e

export PATH="$HOME/.cargo/bin:$PATH"
export RUST_BACKTRACE=1

REPO="$HOME/fast2brad_holo/Fast2bRAD-M"
T2T="/lustre1/g/aos_shihuang/databases/kraken2/kraken16/genomes/GCF_009914755.1_T2T-CHM13v2.0_genomic.fna.gz"
OUTDIR="/home/shihuang/fast2brad_holo/results"

mkdir -p "$OUTDIR"

cd "$REPO"
echo "Building release..."
cargo build --release

echo "Running digests..."
for enzyme in BcgI BsaXI AlfI; do
    echo "=== $enzyme ==="
    target/release/f2brad-host digest \
        -i "$T2T" \
        -s "$enzyme" \
        -o "$OUTDIR/$enzyme" \
        -j 16
done

echo "Done. Results in $OUTDIR"
