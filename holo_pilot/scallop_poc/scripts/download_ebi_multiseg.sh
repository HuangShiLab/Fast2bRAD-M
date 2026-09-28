#!/bin/bash
set -uo pipefail
WORKDIR=/lustre1/g/aos_shihuang/holo2bRAD/scallop_poc
URL="http://ftp.sra.ebi.ac.uk/vol1/fastq/SRR202/008/SRR2027758/SRR2027758.fastq.gz"
FINAL="$WORKDIR/reads/fastq/SRR2027758.fastq.gz"
PARTDIR="$WORKDIR/reads/fastq/parts"
LOG="$WORKDIR/results/multiseg_download.log"
JOBS=4
SEGS=32

echo "$(date): Fetching file size..." | tee "$LOG"
TOTAL=$(curl -sI "$URL" | grep -i content-length | awk "{print \$2}" | tr -d "\r")
if [ -z "$TOTAL" ]; then
    echo "$(date): ERROR could not determine file size" | tee -a "$LOG"
    exit 1
fi
echo "$(date): Total size = $TOTAL bytes" | tee -a "$LOG"

mkdir -p "$PARTDIR"
SEG_SIZE=$((TOTAL / SEGS))

download_seg() {
    local i=$1
    local START=$((i * SEG_SIZE))
    local END=$((START + SEG_SIZE - 1))
    if [ "$i" -eq "$((SEGS - 1))" ]; then
        END=""
    fi
    local PART="$PARTDIR/seg_$i"
    EXPECTED=$(( (i == SEGS - 1) ? TOTAL - START : SEG_SIZE ))

    if [ -f "$PART" ]; then
        PSIZE=$(stat -c%s "$PART" 2>/dev/null || echo 0)
        if [ "$PSIZE" -eq "$EXPECTED" ]; then
            echo "$(date): seg_$i already complete ($PSIZE bytes)" | tee -a "$LOG"
            return 0
        fi
    fi

    echo "$(date): Downloading seg_$i ($START-$END, expected $EXPECTED bytes)" | tee -a "$LOG"
    if [ -z "$END" ]; then
        curl -s --retry 10 --retry-delay 5 -C - -o "$PART" -r "$START-" "$URL"
    else
        curl -s --retry 10 --retry-delay 5 -C - -o "$PART" -r "$START-$END" "$URL"
    fi
    PSIZE=$(stat -c%s "$PART" 2>/dev/null || echo 0)
    echo "$(date): seg_$i finished, size $PSIZE (expected $EXPECTED)" | tee -a "$LOG"
}

export -f download_seg
export URL PARTDIR LOG SEG_SIZE SEGS

echo "$(date): Starting $SEGS segments with up to $JOBS concurrent downloads" | tee -a "$LOG"
for i in $(seq 0 $((SEGS - 1))); do
    download_seg "$i" &
    while [ "$(jobs -r | wc -l)" -ge "$JOBS" ]; do
        wait -n 2>/dev/null || sleep 1
    done
done
wait
echo "$(date): All downloads finished." | tee -a "$LOG"

# verify
MISS=0
for i in $(seq 0 $((SEGS - 1))); do
    EXPECTED=$(( (i == SEGS - 1) ? TOTAL - i * SEG_SIZE : SEG_SIZE ))
    PSIZE=$(stat -c%s "$PARTDIR/seg_$i" 2>/dev/null || echo 0)
    if [ "$PSIZE" -ne "$EXPECTED" ]; then
        echo "$(date): ERROR seg_$i size mismatch ($PSIZE vs $EXPECTED)" | tee -a "$LOG"
        MISS=$((MISS+1))
    fi
done
if [ "$MISS" -ne 0 ]; then
    echo "$(date): $MISS segments incomplete; aborting." | tee -a "$LOG"
    exit 1
fi

echo "$(date): Concatenating segments..." | tee -a "$LOG"
for i in $(seq 0 $((SEGS - 1))); do
    cat "$PARTDIR/seg_$i" >> "$FINAL"
done
echo "$(date): Final file size: $(stat -c%s "$FINAL") bytes" | tee -a "$LOG"

echo "$(date): Decompressing FASTQ..." | tee -a "$LOG"
cd "$WORKDIR/reads/fastq"
pigz -dk "$FINAL" | tee -a "$LOG"
ls -lh | tee -a "$LOG"

echo "$(date): Splitting pooled FASTQ by sample prefix..." | tee -a "$LOG"
python3 - << "PYEOF" | tee -a "$LOG"
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

echo "$(date): Selecting first 3 samples for proof-of-concept..." | tee -a "$LOG"
head -4 "$WORKDIR/reads/by_sample/samples.tsv" > "$WORKDIR/reads/samples.tsv"

echo "$(date): Submitting f2brad-holo classify job..." | tee -a "$LOG"
REPO=/lustre1/g/aos_shihuang/Fast2bRAD-M
HOST_DB=$WORKDIR/db/scallop_BsaXI_host_db.tsv
MICROBE_DB=$WORKDIR/db/BsaXI_example_microbe_db/BsaXI.species.quant.iibdb
MICROBE_DIR=$WORKDIR/db/BsaXI_example_microbe_db
SAMPLES=$WORKDIR/reads/samples.tsv
OUT=$WORKDIR/results

mkdir -p "$OUT"
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
echo "$(date): classify job submitted." | tee -a "$LOG"
