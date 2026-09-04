# Scallop (Wang et al. 2017) holo-2bRAD proof-of-concept

This pilot reproduces the holo-2bRAD concept on *Mizuhopecten yessoensis*
using public BsaXI 2b-RAD data (SRR2027758) and the GCF_002113885.1 reference.

- `scripts/process_after_prefetch.sh` — after SRR2027758.sra is downloaded,
  convert to FASTQ, split by sample prefix, pick the first 3 samples, and
  submit a SLURM `f2brad-holo classify` job.
- `scripts/run_classify.sh` — standalone SLURM job that runs the Rust
  `f2brad-holo classify` binary on `reads/samples.tsv`.
- `scripts/example_microbe_genomes.list` — the 14 example genomes used for the
  minimal BsaXI microbial validation DB.

Data paths on HPC:
- reads: `/lustre1/g/aos_shihuang/holo2bRAD/scallop_poc/reads`
- db:    `/lustre1/g/aos_shihuang/holo2bRAD/scallop_poc/db`
- ref:   `/lustre1/g/aos_shihuang/holo2bRAD/scallop_poc/ref`
- results: `/lustre1/g/aos_shihuang/holo2bRAD/scallop_poc/results`

## Fallback: EBI FASTQ download
If the NCBI SRA prefetch stalls, `scripts/download_ebi_then_process.sh` fetches the
same run directly from ENA (`SRR2027758.fastq.gz`), decompresses with `pigz`, splits
by sample prefix, and submits the SLURM `f2brad-holo classify` job.
