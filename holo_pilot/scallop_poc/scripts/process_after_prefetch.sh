#!/bin/bash
set -euo pipefail

WORKDIR=/lustre1/g/aos_shihuang/holo2bRAD/scallop_poc
REPO=/lustre1/g/aos_shihuang/Fast2bRAD-M
SRA=$WORKDIR/reads/sra/SRR2027758/SRR2027758.sra
OUT=$WORKDIR/results

echo "$(date): Waiting for prefetch to complete..."
while [ ! -f "$SRA" ]; do
    if [ -f "$SRA.tmp" ]; then
        size=$(du -h "$SRA.tmp" 2>/dev/null | cut -f1)
        echo "$(date): Downloading... $size"
    fi
    sleep 120
done
echo "$(date): SRA file ready: $SRA"

echo "Converting to FASTQ with fasterq-dump..."
mkdir -p "$WORKDIR/reads/fastq"
cd "$WORKDIR/reads/fastq"
$REPO/../tools/sratoolkit.3.0.0-centos_linux64/bin/fasterq-dump --split-files -O . "$SRA" 2>&1 | tail -20
ls -lh

echo "Splitting pooled FASTQ by sample prefix..."
python3 - << "PYEOF"
import gzip, re, os
from collections import defaultdict

fastq = "/lustre1/g/aos_shihuang/holo2bRAD/scallop_poc/reads/fastq/SRR2027758.fastq"
outdir = "/lustre1/g/aos_shihuang/holo2bRAD/scallop_poc/reads/by_sample"
os.makedirs(outdir, exist_ok=True)

handles = {}
meta = []

def get_handle(sample):
    if sample not in handles:
        path = os.path.join(outdir, f"{sample}.fastq.gz")
        handles[sample] = gzip.open(path, "wt")
        meta.append((sample, path))
    return handles[sample]

with open(fastq) as fh:
    i = 0
    for line in fh:
        if i % 4 == 0:
            # @F13-2_L8_I001:... or similar
            m = re.match(r"@([^:_\s]+)", line)
            sample = m.group(1) if m else "unknown"
            get_handle(sample).write(line)
        else:
            get_handle(sample).write(line)
        i += 1
        if i % 4000000 == 0:
            print(f"Processed {i//4} reads")

for h in handles.values():
    h.close()

with open(os.path.join(outdir, "samples.tsv"), "w") as out:
    out.write("sample_id\tr1\n")
    for sample, path in sorted(meta):
        out.write(f"{sample}\t{path}\n")

print(f"Done. {len(handles)} samples written to {outdir}")
PYEOF

echo "Selecting first 3 samples for proof-of-concept..."
head -4 "$WORKDIR/reads/by_sample/samples.tsv" > "$WORKDIR/reads/samples.tsv"

echo "Submitting f2brad-holo classify job..."
cat > "$OUT/classify_job.sh" << "EOFCLASS"
#!/bin/bash
#SBATCH --job-name=scallop_classify
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=04:00:00
#SBATCH --output=$OUT/classify_%j.log
set -euo pipefail

WORKDIR=/lustre1/g/aos_shihuang/holo2bRAD/scallop_poc
REPO=/lustre1/g/aos_shihuang/Fast2bRAD-M
HOST_DB=$WORKDIR/db/scallop_BsaXI_host_db.tsv
MICROBE_DB=$WORKDIR/db/BsaXI_example_microbe_db/BsaXI.species.quant.iibdb
MICROBE_DIR=$WORKDIR/db/BsaXI_example_microbe_db
SAMPLES=$WORKDIR/reads/samples.tsv
OUTDIR=$WORKDIR/results/holo

mkdir -p "$OUTDIR"

$REPO/target/release/f2brad-holo classify \\
    --host-db "$HOST_DB" \\
    --microbe-db "$MICROBE_DB" \\
    --microbe-db-dir "$MICROBE_DIR" \\
    --site BsaXI \\
    --sample-list "$SAMPLES" \\
    --output "$OUTDIR" \\
    --taxonomy species \\
    --host-max-mismatch 2 \\
    --min-depth 4 \\
    -j 16

echo "Done. Results in $OUTDIR"
EOFCLASS

sbatch "$OUT/classify_job.sh"
echo "$(date): classify job submitted."
