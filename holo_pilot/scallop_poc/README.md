# Scallop (Wang et al. 2017) holo-2bRAD proof-of-concept

This pilot reproduces the holo-2bRAD concept on *Mizuhopecten yessoensis*
using public BsaXI 2b-RAD data (SRR2027758) and the GCF_002113885.1 reference.

## Scripts
- `scripts/process_after_prefetch.sh` — after SRR2027758.sra is downloaded,
  convert to FASTQ, split by sample prefix, pick the first 3 samples, and
  submit a SLURM `f2brad-holo classify` job.
- `scripts/run_classify.sh` — standalone SLURM job that runs the Rust
  `f2brad-holo classify` binary on `reads/samples.tsv`.
- `scripts/download_ebi_then_process.sh` — single-stream EBI FASTQ fallback.
- `scripts/download_ebi_multiseg.sh` — multi-segment HTTP range downloader
  with resume/retry (used when NCBI prefetch stalled).
- `scripts/example_microbe_genomes.list` — the 14 example genomes used for the
  minimal BsaXI microbial validation DB.

## Important data note
SRR2027758 from ENA is the pooled 2b-RAD library with SRA-style read names
(`@SRR2027758.N ...`). The original per-sample barcodes are not present in the
FASTQ headers, so the pipeline was run on the whole pool (subsampled to 200 M
reads for the first proof-of-concept). Host genotypes therefore represent a
population consensus, not an individual.

## PoC results (200 M read subsample)
See `results/scallop_poc_summary.md` on the HPC and the synced copy in
`Fast2bRAD-M-paper/holo2bRAD/results/scallop_poc/`.

Highlights:
- Host-tag mapping rate: ~36.6% of input fragments (73.2 M / 200 M)
- Host fraction of classified fragments: ~1.0
- Variant SNP sites: ~360 k
- Ti/Tv: 0.87 (low because the library is pooled)

## HPC paths
- reads: `/lustre1/g/aos_shihuang/holo2bRAD/scallop_poc/reads`
- db:    `/lustre1/g/aos_shihuang/holo2bRAD/scallop_poc/db`
- ref:   `/lustre1/g/aos_shihuang/holo2bRAD/scallop_poc/ref`
- results: `/lustre1/g/aos_shihuang/holo2bRAD/scallop_poc/results`
