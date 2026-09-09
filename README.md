# Fast2bRAD-M

**Fast2bRAD-M** is a high-performance Rust reimplementation of the [2bRAD-M](https://github.com/HuangShiLab/2bRAD-M) microbiome profiling pipeline. It delivers the same analytical results as the original Perl/Shell pipeline while achieving dramatically higher throughput through native parallelism and optimized I/O.

---

## Table of Contents

- [Features](#features)
- [Installation](#installation)
- [Quick Start](#quick-start)
- [Pipeline Overview](#pipeline-overview)
- [Subcommands](#subcommands)
  - [fast2bRAD-M core](#fast2brad-m-core)
    - [extract](#extract)
    - [build-qual-db](#build-qual-db)
    - [build-quan-db](#build-quan-db)
    - [dedup-db](#dedup-db)
    - [quantify](#quantify)
    - [find-genome](#find-genome)
    - [merge](#merge)
    - [predict](#predict)
    - [inspect](#inspect)
    - [pipeline](#pipeline)
  - [f2brad-host](#f2brad-host)
    - [digest](#digest)
    - [cross](#cross)
    - [build-db](#build-db)
    - [genotype](#genotype)
  - [f2brad-holo](#f2brad-holo)
    - [classify](#classify)
- [File Formats](#file-formats)
- [Supported Enzymes](#supported-enzymes)
- [Output Directory Structure](#output-directory-structure)
- [License](#license)

---

## Features

- **High Performance** — Rust implementation with Rayon multi-core parallelism; batch-digesting 15 reference genomes in < 0.12 s
- **Full Enzyme Support** — All 16 Type IIB restriction enzymes (BcgI, CspCI, AloI, BsaXI, BaeI, CjeI, PpiI, PsrI, BplI, FalI, Bsp24I, HaeIV, CjePI, Hin4I, AlfI, BslFI)
- **All Input Types** — Reference genomes, Shotgun metagenomic reads (SE/PE), and single 2bRAD tags
- **Built-in QC** — N-ratio, minimum quality score, and minimum quality-percent filtering
- **Functional Prediction** — Matrix-multiplication-based functional abundance profiling (KO, KEGG, etc.)
- **Host Genotyping** — `f2brad-host` builds a host tag database and calls genotypes from 2bRAD reads
- **Holo-2bRAD Integration** — `f2brad-holo` performs one-pass joint host genotyping + microbial profiling with microbial cross-assignment masking
- **Viral Profiling (VIP2B)** — Convert the VIP2B UHGV viral database and profile viruses with the same 8-enzyme strategy
- **Database Inspection** — `inspect` reports format, tag counts and example records from binary `.iibdb`/`.iibsp` files
- **Resume Support** — `.done` marker files allow interrupted runs to be resumed without re-computation
- **One-Command Pipeline** — The `pipeline` subcommand chains all steps automatically

---

## Installation

### Option 1 — Conda (Recommended)

```bash
conda env create -f fast2brad_m_conda.yaml -n fast2brad
conda activate fast2brad
cargo build --release
```

### Option 2 — Direct Compilation

Prerequisites: [Rust toolchain](https://rustup.rs/) ≥ 1.70

```bash
git clone https://github.com/HuangShiLab/Fast2bRAD-M.git
cd Fast2bRAD-M
cargo build --release
# Binary: target/release/fast2bRAD-M
```

> **Note**: Paired-end PEAR merging (optional) requires PEAR to be installed separately:
> ```bash
> conda install -c bioconda pear
> ```

---

## Quick Start

```bash
# One-command full pipeline (database construction + sample profiling)
fast2bRAD-M pipeline \
  --mode full \
  --samples samples.tsv \
  --genome-list genome_list.tsv \
  --taxonomy taxonomy.tsv \
  --site BcgI \
  --level species \
  --outdir results/ \
  --prefix my_run \
  --threads 16 \
  --resume yes

# With functional prediction
fast2bRAD-M pipeline \
  --mode full \
  --samples samples.tsv \
  --taxonomy taxonomy.tsv \
  --site BcgI \
  --level species \
  --outdir results/ \
  --prefix my_run \
  --threads 16 \
  --ko-mapping ko_mapping.tsv
```

---

## Pipeline Overview

The full analysis pipeline runs in five main stages:

```
Raw reads (FASTQ)
      │
      ▼
[1] extract          →  01_extract/{sample}.BcgI.iibsp
      │
      ▼
[2] build-qual-db    →  02_db_qual/  (qualitative database, shared)
      │
      ▼
[3] quantify (qual)  →  qualitative/{sample}/  (qualitative screen)
      │
      ▼
[4] find-genome      →  quantitative_sdb/{sample}/sdb.list
      │
      ▼
[5] build-quan-db  } →  02_db_quan/{sample}/  (per-sample quantitative DB)
    quantify (quan) } →  04_quantify/{sample}/
      │
      ▼
[6] merge            →  05_merge/{prefix}.all.xls
      │
      ▼ (optional, requires --ko-mapping)
[7] predict          →  05_merge/{prefix}.func.xls
```

`dedup-db` provides an alternative way to obtain a quantitative database: it converts a qualitative `.iibdb` into a quantitative one by dropping any tag that maps to more than one GCF. This single-GCF database is the format consumed by `f2brad-holo classify`.

---

## Subcommands

The project provides three command-line binaries:

| Binary | Purpose |
|--------|---------|
| `fast2bRAD-M` | Core microbiome profiling pipeline (extract → build-db → quantify → merge → predict) |
| `f2brad-host` | Host 2bRAD analysis: in-silico digest, microbial cross-assignment masking, host DB construction, and genotyping |
| `f2brad-holo` | One-pass holo-2bRAD driver: joint host genotyping + microbial profiling |

---

## fast2bRAD-M core

### `extract`

Digest input sequences with a Type IIB restriction enzyme and extract 2bRAD tags.

```bash
fast2bRAD-M extract \
  --genome-list sample_list.tsv \  # batch mode
  -t 2 \                           # input type (1=reference, 2=shotgun, 3=single tag)
  -s BcgI \                        # enzyme name or ID (1–16)
  --od output_dir/ \
  --op sample_prefix \
  -j 8 \                           # threads
  --qc yes \                       # quality control
  --qc-scope auto \                # QC on the whole read or on the tag window
  -n 0.08 \                        # max N ratio
  -q 30 \                          # min quality score
  -p 80                            # min quality percent
```

`--qc-scope` decides which bases the `-n/-q/-p` filters look at: `read` drops a
read whose low-quality or N bases lie anywhere in it, `tag` judges only the tag
window so distant bad bases in a long shotgun read no longer discard an
otherwise clean tag. `auto` (the default) keeps each input type's historical
behaviour — `read` for `-t 2`, `tag` for `-t 3`.

**Input types**:

| Type | Description |
|------|-------------|
| 1    | Reference genome FASTA (for database construction) |
| 2    | Shotgun metagenome reads (SE or PE; PE can use PEAR merging) |
| 3    | Single 2bRAD tag reads |

**Output**:
- `{prefix}.{enzyme}.iibsp` — Binary tag file for sample reads (Types 2 & 3)
- `{prefix}.{enzyme}.iibdb` — Binary tag file for reference genomes (Type 1)
- `{prefix}.{enzyme}.stat.tsv` — Digest statistics. For `-t 1` the columns are
  `sample/enzyme/contigs/total_bases/tag_count/tags_per_mb`; for `-t 2`/`-t 3`
  they are read-based counts

**Paired-end input** (two files given to `-i`, or a third column in the `-l`
list): both mates are always processed. Without PEAR they are read in lockstep
and a tag seen on both mates of a pair is counted once, since overlapping mates
describe the same physical fragment.

**Batch mode failures**: a sample that fails to process has its partial output
removed and the command exits non-zero, listing the failed samples — an empty
or truncated database is never left behind for downstream steps to misread as a
tag-free sample.

**Paired-end with PEAR merging** (optional, Type 2 only):
```bash
fast2bRAD-M extract \
  -i sample_R1.fq.gz sample_R2.fq.gz \
  -t 2 -s BcgI \
  --od output/ --op sample1 \
  --use-pear yes --pe pear --pc 4
```

---

### `build-qual-db`

Build a qualitative (classification-specificity) database from reference genomes.

```bash
fast2bRAD-M build-qual-db \
  -l genome_list.tsv \   # genome list (2-column: genome_id + fasta_path)
  --taxonomy taxonomy.tsv \ # taxonomy file (genome_id + taxonomy columns)
  -s BcgI \              # enzyme
  -t species \           # taxonomy level(s); comma-separated or "all"
  -o db_qual/ \
  --pre-digested-dir pre_digested/ \  # optional: pre-digested .iibdb files
  -r yes \               # remove redundant tags
  -j 8
```

**Genome list format** (2-column, tab-separated):
```
GCF_000007445.1  /path/to/genome.fna.gz
GCF_000007445.2  /path/to/another_genome.fna.gz
```

**Taxonomy file format** (tab-separated, 9 columns):
```
GCF_000007445.1  Bacteria  Proteobacteria  Gammaproteobacteria  Enterobacterales  Enterobacteriaceae  Escherichia  Escherichia_coli  str.K-12
```
Or GTDB format (second column = `d__Bacteria;p__Proteobacteria;...`).

Backward compatibility: If `--taxonomy` is not provided, `--list` can also be a single file with both genome paths and taxonomy (original format).

**Output** (per taxonomy level):
- `{enzyme}.enzyme.iibdb` — All tags from all genomes (intermediate)
- `{enzyme}.{level}.iibdb` — Taxon-unique tags only
- `abfh_classify_with_speciename.txt.gz` — GCF-to-taxonomy mapping

---

### `build-quan-db`

Build a quantitative (per-sample) database that retains only unique tags.

```bash
fast2bRAD-M build-quan-db \
  -l sdb.list \           # genome list for this sample (from find-genome)
  -s BcgI \
  -t species \
  -o sample_db/ \
  -e qual_db/BcgI.enzyme.iibdb \  # reuse the enzyme file from qual DB
  -j 4
```

**Output**:
- `BcgI.species.iibdb` — Unique tags for quantitative profiling
- `abfh_classify_with_speciename.txt.gz` — Taxonomy mapping

---

### `quantify`

Calculate per-taxon relative abundance for one or more samples.

```bash
fast2bRAD-M quantify \
  -l sample_list.tsv \   # sample_name<TAB>path_to.iibsp
  -d database_dir/ \     # directory with BcgI.species.iibdb + taxonomy mapping file
  -t species \
  -s BcgI \
  -o quantify_out/ \
  -g 5.0 \               # G-score threshold (species with G < threshold excluded)
  -v yes \               # verbose: output per-tag detail files
  -j 8
```

**G-score** = `sqrt(sequenced_tag_num × sequenced_reads_num)` — a combined measure of breadth and depth of coverage.

**Output** per sample (inside `output_dir/{sample}/`):
- `{sample}.{enzyme}.xls` — Per-taxon abundance table with statistics
- `{sample}.{enzyme}.GCF_detected.xls` — Per-genome detection details

**Abundance table columns**:
```
Kingdom  Phylum  Class  Order  Family  Genus  Species
Theoretical_Tag_Num  Sequenced_Tag_Num  Percent
Sequenced_Reads_Num  Reads/Theoretical  Reads/Sequenced
Sequenced_Tag_Num(depth>1)  G_Score
```

---

### `find-genome`

Filter reference genomes for quantitative analysis based on qualitative results.
This step converts broad qualitative detections into a per-sample genome list.

```bash
fast2bRAD-M find-genome \
  -l samples.tsv \
  -d qual_db/ \
  -o quantitative_sdb/ \
  --qual-dir qualitative/ \
  --gscore 5 \     # G-score threshold for qualitative detection
  --gcf 1 \        # minimum detected tags per GCF
  -j 8
```

**Output** per sample: `quantitative_sdb/{sample}/sdb.list` — tab-separated genome records that pass thresholds.

---

### `merge`

Merge per-sample quantitative results into a combined abundance table.

```bash
fast2bRAD-M merge \
  -l merge_list.tsv \   # sample_name<TAB>path_to_{sample}.{enzyme}.xls
  -o merge_out/ \
  -p Abundance_Stat \   # output file prefix
  --mock mock1,mock2 \  # comma-separated mock sample names (filtered out)
  --control ctrl1       # comma-separated negative control names
```

**Output**:
- `{prefix}.all.xls` — Merged relative abundance matrix (all samples)
- `{prefix}.filtered.xls` — Same, with mock/control samples and contamination taxa removed

**Merge table format**:
```
Kingdom  Phylum  Class  Order  Family  Genus  Species  sample1  sample2  ...
Bacteria  Proteobacteria  ...  Escherichia_coli  0.3413  0.2841  ...
```
Values are relative abundances normalized to sum to 1.0 per sample.

---

### `predict`

Predict functional abundance by multiplying the species abundance matrix with a species-to-function mapping matrix.

**Formula**: `Functional_abundance = t(Species_abundance) × Mapping_matrix`

```bash
fast2bRAD-M predict \
  -a 05_merge/Abundance_Stat.all.xls \   # merged species abundance table
  -m ko_mapping.tsv \                    # species-to-KO mapping matrix
  -o 05_merge/ \
  -p Abundance_Stat
```

**Mapping matrix format** (TSV):
```
#Species         KO00001  KO00002  KO00003  ...
Escherichia_coli   5        0        3       ...
Cutibacterium_acnes 2       8        0       ...
```
- First column: species name (must match the Species column in the abundance table)
- Remaining columns: KO/functional IDs; values = gene copy counts

**Output**:
- `{prefix}.func.xls` — Functional abundance table, per-sample normalized (each sample sums to 1.0)

```
#Function  sample1     sample2     ...
KO00001    0.12345678  0.09876543  ...
KO00002    0.00000000  0.04321098  ...
```

---

### `dedup-db`

Convert a qualitative compact database (`.iibdb`) into a quantitative database by keeping only tags that map to a single GCF. This is the database format expected by `f2brad-holo classify` for microbial profiling.

```bash
fast2bRAD-M dedup-db \
  -i qual_db/BcgI.species.iibdb \
  -o quan_db/BcgI.species.quant.iibdb
```

**Parameters**:
| Parameter | Required | Description |
|-----------|----------|-------------|
| `-i` / `--input` | Yes | Input qualitative `.iibdb` |
| `-o` / `--output` | Yes | Output quantitative `.iibdb` |

**Output**:
- A quantitative `.iibdb` where every tag hash is unique to one reference genome

---

### `inspect`

Inspect `.iibdb` / `.iibsp` binary files: show format, tag counts and example records.

```bash
fast2bRAD-M inspect qual_db/BcgI.species.iibdb

# Full scan with per-genome counts
fast2bRAD-M inspect -f -t 20 qual_db/BcgI.species.iibdb
```

**Parameters**:
| Parameter | Default | Description |
|-----------|---------|-------------|
| `-f` / `--full` | off | Count every tag (and per-genome tag counts). Slow on large DBs |
| `--distinct` | off | With `--full`, also count distinct tag hashes (memory-heavy) |
| `-n` / `--records` | 5 | Number of example records to print per file |
| `-t` / `--top` | 10 | Number of top genomes to list |
| `-o` / `--output` | stdout | Write report to file |

**Output**:
- Human-readable report: file format, record count, example records, and (with `--full`) per-genome tag counts

---

### `pipeline`

One-command orchestrator that chains all steps automatically.

#### Run Modes

| Mode | Description |
|------|-------------|
| `full` | Build database + profile all samples |
| `db-only` | Build qualitative database only |
| `sample-only` | Profile samples using an existing database |

#### Full Pipeline

```bash
fast2bRAD-M pipeline \
  --mode full \
  --samples samples.tsv \
  --genome-list genome_list.tsv \
  --taxonomy taxonomy.tsv \
  --site BcgI \
  --level species \
  --outdir results/ \
  --prefix run1 \
  --threads 16 \
  --gscore 5 \
  --gcf 1 \
  --resume yes
```

#### Database Build Only

```bash
fast2bRAD-M pipeline \
  --mode db-only \
  --genome-list genome_list.tsv \
  --taxonomy taxonomy.tsv \
  --pre-digested-dir pre_digested/ \
  --site BcgI \
  --level species \
  --outdir db/ \
  --threads 16
```

#### Sample-Only (Use Existing Database)

```bash
fast2bRAD-M pipeline \
  --mode sample-only \
  --samples samples.tsv \
  --database db/ \
  --site BcgI \
  --level species \
  --outdir results/ \
  --prefix run1 \
  --threads 16 \
  --resume yes
```

#### With Functional Prediction

```bash
fast2bRAD-M pipeline \
  --mode sample-only \
  --samples samples.tsv \
  --database db/ \
  --site BcgI \
  --outdir results/ \
  --prefix run1 \
  --threads 16 \
  --ko-mapping ko_mapping.tsv   # triggers automatic predict step after merge
```

#### All Pipeline Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--mode` | `full` | Run mode: `full`, `db-only`, `sample-only` |
| `--samples` / `-l` | — | Sample list TSV (required for `full`/`sample-only`) |
| `--genome-list` | — | Reference genome list (for `db-only` / database building) |
| `--taxonomy` | — | Taxonomy/classify file (TSV or GTDB format) |
| `--database` | — | Pre-built database directory (for `sample-only`) |
| `--pre-digested-dir` | — | Directory with pre-digested `.iibdb` files |
| `--site` / `-s` | — | Enzyme name (`BcgI`) or ID (`1`–`16`) |
| `--level` / `-t` | `species` | Taxonomy level for profiling |
| `--outdir` | — | Output directory |
| `--prefix` | `Abundance_Stat` | Prefix for output files |
| `--threads` / `-j` | auto | Global thread count |
| `--gscore` | `5.0` | G-score threshold for find-genome |
| `--gcf` | `1` | Min detected tags per GCF in find-genome |
| `--resume` | `no` | Skip steps that already have `.done` markers (`yes`/`no`) |
| `--qc` | `yes` | Quality control for extract |
| `--max-n` | `0.08` | Max N-base ratio |
| `--min-qual` | `30` | Min base quality score |
| `--min-qual-percent` | `80` | Min percent of bases passing quality |
| `--qual-base` | `33` | Quality score encoding base |
| `--use-pear` | `no` | Enable PEAR merging for paired-end reads |
| `--pear-bin` | `pear` | Path to PEAR executable |
| `--pc` | `1` | Threads per PEAR process |
| `--mock` | — | Comma-separated mock sample names (for merge filtering) |
| `--control` | — | Comma-separated negative control names (for merge filtering) |
| `--ko-mapping` | — | Species-to-function mapping matrix; triggers `predict` step after merge |

---

## `f2brad-host`

Host-side utilities for holo-2bRAD analysis. The typical workflow is:

1. `digest` the host reference genome (e.g. T2T-CHM13v2.0) to obtain per-locus tags.
2. `cross` compare those human tags against a microbial genome database to identify tags that could be mis-assigned to microbes.
3. `build-db` create a masked host tag database, optionally removing cross-assignable tags.
4. `genotype` a 2bRAD sample against the host database.

### `digest`

In-silico digest a reference genome and report tag-level statistics.

```bash
f2brad-host digest \
  -i chm13v2.0.fa.gz \
  -s BcgI \
  -o chm13v2.0_BcgI_digest/ \
  -j 8
```

**Output**:
- `sites.tsv` — one row per tag locus with sequence, canonical sequence, hash, GC fraction, CpG count, and uniqueness flag
- `stat.tsv` — summary statistics

### `cross`

Cross-assignment collision analysis: scan human tags against microbial genomes to find tags that match microbial sequences within a Hamming-distance threshold. These tags should be masked from the host genotype database when analyzing human microbiome samples.

```bash
f2brad-host cross \
  -t chm13v2.0_BcgI_digest/sites.tsv \
  -l microbial_genome_list.tsv \
  -s BcgI \
  -o chm13v2.0_BcgI_cross/ \
  --max-mismatch 2 \
  -j 16
```

**Output**:
- `collisions.tsv` — per-human-tag collision report
- `mask.list` — one canonical hash per line, ready for `build-db --human-mask`

### `build-db`

Build a host tag database from a `digest` sites file and an optional cross-assignment mask.

```bash
f2brad-host build-db \
  -t chm13v2.0_BcgI_digest/sites.tsv \
  -m chm13v2.0_BcgI_cross/mask.list \
  -s BcgI \
  -o chm13v2.0_BcgI.host_db.tsv
```

**Output**:
- A TSV host tag database with columns `contig`, `pos`, `strand`, `seq`, `canonical`, `hash`, `gc_frac`, `cpg_count`, `cpg_island`, `unique`

### `genotype`

Genotype a 2bRAD sample against a host tag database.

```bash
f2brad-host genotype \
  -d chm13v2.0_BcgI.host_db.tsv \
  -1 sample_R1.fq.gz \
  -2 sample_R2.fq.gz \
  -s BcgI \
  -o sample_genotype/ \
  --max-mismatch 2 \
  --min-depth 4 \
  -j 8
```

**Output**:
- `genotypes.vcf` — per-locus diploid genotype calls (GT/DP/AD/PL) in VCF 4.2 format
- `dosages.bimbam` — mean dosages for downstream SNP-based analyses

---

## `f2brad-holo`

One-pass holo-2bRAD driver that jointly profiles host genotypes and microbial composition from the same 2bRAD library.

### `classify`

Run host genotyping and microbial profiling in a single pass. This is useful for human microbiome samples where the same sequencing reads contain both host and microbial 2bRAD tags.

```bash
f2brad-holo classify \
  -d chm13v2.0_BcgI.host_db.tsv \
  -m microbial_db/BcgI.species.quant.iibdb \
  --microbe-db-dir microbial_db/ \
  --microbe-mask chm13v2.0_BcgI_cross/mask.list \
  -1 sample_R1.fq.gz \
  -2 sample_R2.fq.gz \
  -s BcgI \
  -o holo_results/sample1/ \
  --sample-name sample1 \
  --exclude-human \
  -j 8
```

Batch mode (process many samples in parallel):

```bash
f2brad-holo classify \
  -d chm13v2.0_BcgI.host_db.tsv \
  -m microbial_db/BcgI.species.quant.iibdb \
  --microbe-db-dir microbial_db/ \
  --microbe-mask chm13v2.0_BcgI_cross/mask.list \
  -l samples.tsv \
  -s BcgI \
  -o holo_results/ \
  --exclude-human \
  -j 16
```

**Parameters**:
| Parameter | Required | Description |
|-----------|----------|-------------|
| `-d` / `--host-db` | Yes | Host tag database TSV from `f2brad-host build-db` |
| `-m` / `--microbe-db` | Yes | Microbial quantitative `.iibdb` (single-GCF, from `dedup-db`) |
| `--microbe-db-dir` | No | Directory with `{enzyme}.{level}.iibdb` + `abfh_classify_with_speciename.txt.gz`; enables `species_counts.tsv` |
| `--microbe-mask` | No | Canonical-tag mask of human-to-microbe cross-assignable tags |
| `-1` / `--r1` | Yes* | Read 1 FASTQ (gzip ok). *Ignored when `-l` is used |
| `-2` / `--r2` | No | Read 2 FASTQ (gzip ok). *Ignored when `-l` is used |
| `-l` / `--sample-list` | No | Batch sample list: `sample_name<TAB>r1_path[<TAB>r2_path]` |
| `-s` / `--site` | Yes | Enzyme name or ID (1–16) |
| `-o` / `--output` | Yes | Output directory |
| `--sample-name` | No | Sample name for single-sample mode [default: `sample`] |
| `--host-max-mismatch` | No | Max Hamming distance for host tag matching [default: 2] |
| `-q` / `--min-qual` | No | Min Phred quality for genotype pileup [default: 20] |
| `--min-depth` | No | Min per-locus depth to emit a genotype [default: 4] |
| `-t` / `--taxonomy` | No | Taxonomy level for microbial counts [default: `species`] |
| `--output-iibsp` | No | Write `sample.iibsp.gz` for downstream `fast2bRAD-M quantify` |
| `--exclude-human` | No | Drop GCFs whose species is `human` from the microbial DB |
| `-j` / `--threads` | No | Threads [default: 4] |

**Output**:
- `genotypes.vcf` — host genotype calls
- `species_counts.tsv` — microbial taxon counts (when `--microbe-db-dir` is provided)
- `holo_classify.tsv` — read-classification summary (host/microbe/ambiguous fractions)
- `sample.iibsp.gz` — optional sample tag stream for downstream `fast2bRAD-M quantify` (with `--output-iibsp`)

---

## Viral profiling

Fast2bRAD-M ships with helper scripts that reuse 8-enzyme viral databases inside the existing `fast2bRAD-M quantify` framework. Two databases are currently supported:

- **[VIP2B](https://github.com/sunzhengCDNM/VIP2B)** — gut/phage database based on UHGV
- **[HOVD](https://hovd.org)** — Human Oral Virome Database (OPD + OED)

Both use the same 8-enzyme strategy (AlfI, BcgI, BslFI, CjeI, CjePI, FalI, HaeIV, Hin4I) and the same `quantify_vip2b.py` wrapper.

### VIP2B database

#### 1. Obtain the VIP2B database

Download the three files below from the VIP2B Zenodo record (e.g. `https://zenodo.org/records/18944630`):

- `8Enzyme.Species.uniq.marisa`
- `abfh_classify_with_speciename.txt.gz`
- `metadata.tsv.gz`

#### 2. Convert to fast2bRAD-M format

```bash
python tools/convert_vip2b_db.py \
  -m 8Enzyme.Species.uniq.marisa \
  -c abfh_classify_with_speciename.txt.gz \
  -o vip2b_db/ \
  -l species \
  -s BcgI
```

This produces:
- `vip2b_db/BcgI.species.iibdb` (placeholder enzyme name for the combined 8-enzyme DB)
- `vip2b_db/abfh_classify_with_speciename.txt.gz`
- `vip2b_db/BcgI.species.iibdb.stats.txt`

> **Note on the `-s` placeholder**: VIP2B's `8Enzyme.Species.uniq` database already contains tags from all eight enzymes. `fast2bRAD-M quantify` requires a valid enzyme name purely for output naming, so `convert_vip2b_db.py` defaults to `BcgI`. The actual quantification below extracts tags with all eight enzymes and merges them before profiling.

### HOVD database

#### 1. Obtain HOVD

Download the genome FASTA and annotation table from [https://hovd.org](https://hovd.org):

- `HOVD-geneseqences.fasta`
- `HOVD-annotations.xlsx`

#### 2. Build the 8-enzyme fast2bRAD-M database

```bash
python tools/build_hovd_db.py \
  --fasta HOVD-geneseqences.fasta \
  --annotation HOVD-annotations.xlsx \
  -o hovd_db/ \
  -j 8
```

This produces:
- `hovd_db/BcgI.species.iibdb`
- `hovd_db/abfh_classify_with_speciename.txt.gz`
- `hovd_db/metadata.tsv.gz`
- `hovd_db/BcgI.species.iibdb.stats.txt`

> **Taxonomy note**: HOVD marks many OPD contigs with `uc_*` (unclassified) at species/genus level. `build_hovd_db.py` converts these to `unknown` and uses the HOVD `contig_id` as the `Species` column so that every abundance row remains identifiable and can be annotated with metadata. You can collapse the resulting contig-level profile to family/class/phylum using `annotate_hovd.py` or your own post-processing.

### Profile samples

Prepare a sample list (`samples.tsv`):

```tsv
sample1  /path/sample1_R1.fq.gz  /path/sample1_R2.fq.gz
sample2  /path/sample2_R1.fq.gz
```

Run the 8-enzyme wrapper (replace `hovd_db/` with `vip2b_db/` for VIP2B):

```bash
python tools/quantify_vip2b.py \
  -i samples.tsv \
  -d hovd_db/ \
  -o hovd_results/ \
  -j 16 \
  --qc no \
  --merge-prefix HOVD
```

Outputs:
- `{outdir}/01_extract/` — per-enzyme `.iibsp` files and combined `.VIP2B.iibsp`
- `{outdir}/02_quantify/{sample}/{sample}.BcgI.xls` — per-sample viral abundance
- `{outdir}/03_profiles/{sample}.VIP2B.xls` — convenience copy of each profile
- `{outdir}/{prefix}.all.xls` and `{prefix}.filtered.xls` — merged abundance matrix (if ≥2 samples)

For already-demultiplexed 2bRAD tag reads, use `--input-type 3`. For WGS/shotgun data (default), use `--input-type 2`.

### Annotate with viral metadata

For VIP2B:

```bash
python tools/annotate_vip2b.py \
  -i vip2b_results/VIP2B.all.xls \
  -d vip2b_db/metadata.tsv.gz \
  -o vip2b_results/annotation/
```

This writes:
- `Phenotype.tsv` — lifestyle, jumbo-phage status, viralverify prediction
- `viral_function/Uniref90.tsv` — normalized UniRef90 gene abundance
- `viral_function/cluster.tsv` — normalized info-annotation cluster abundance
- `viral_taxonomy/{kingdom..genus}_abund.tsv` — collapsed viral taxonomy
- `host_taxonomy/{domain..species}_abund.tsv` — collapsed host taxonomy

For HOVD:

```bash
python tools/annotate_hovd.py \
  -i hovd_results/HOVD.all.xls \
  -d hovd_db/metadata.tsv.gz \
  -o hovd_results/annotation/
```

This writes:
- `metadata_stats.tsv` — checkv quality, provirus status, geography, oral site
- `viral_taxonomy/{kingdom..species}_abund.tsv` — collapsed oral viral taxonomy

---

## File Formats

### Sample List (`samples.tsv`)

```tsv
# sample_name  path_to_R1               path_to_R2 (optional for PE)
sample1         /path/sample1_R1.fq.gz  /path/sample1_R2.fq.gz
sample2         /path/sample2_R1.fq.gz
```

### Genome List (`genome_list.tsv`)

Standard format:
```tsv
GCF_000007445.1  Bacteria  Proteobacteria  Gammaproteobacteria  Enterobacterales  Enterobacteriaceae  Escherichia  Escherichia_coli  str.K-12  /path/to/genome.fna.gz
```

GTDB format (auto-detected):
```tsv
GCF_000007445.1  d__Bacteria;p__Proteobacteria;c__Gammaproteobacteria;...
```

### KO Mapping Matrix (`ko_mapping.tsv`)

```tsv
#Species                KO00001  KO00002  KO00003
Escherichia_coli           5        0        3
Cutibacterium_acnes        2        8        0
```

---

## Supported Enzymes

Tag lengths below are derived directly from the `@site` regex patterns in
the original `2bRADExtraction.pl` (shihuang047/2bRAD-M), so extracted tags
match the Perl implementation exactly for every enzyme, not just BcgI.
BaeI, HaeIV, and Hin4I have IUPAC-degenerate positions in their recognition
sequence (e.g. a pyrimidine-only or purine-only base) that cannot be
represented as a fixed literal byte string, so they are matched with an
anchored regex instead of the fixed-byte-offset matcher used for the other
13 enzymes; matching behavior is otherwise identical.

| ID | Name    | Tag Length | Matching mode |
|----|---------|-----------|----------------|
| 1  | CspCI   | 33 bp | fixed-byte |
| 2  | AloI    | 27 bp | fixed-byte |
| 3  | BsaXI   | 27 bp | fixed-byte |
| 4  | BaeI    | 28 bp | regex (degenerate base) |
| **5**  | **BcgI** *(recommended)* | **32 bp** | fixed-byte |
| 6  | CjeI    | 28 bp | fixed-byte |
| 7  | PpiI    | 27 bp | fixed-byte |
| 8  | PsrI    | 27 bp | fixed-byte |
| 9  | BplI    | 27 bp | fixed-byte |
| 10 | FalI    | 27 bp | fixed-byte |
| 11 | Bsp24I  | 27 bp | fixed-byte |
| 12 | HaeIV   | 27 bp | regex (degenerate base) |
| 13 | CjePI   | 27 bp | fixed-byte |
| 14 | Hin4I   | 27 bp | regex (degenerate base) |
| 15 | AlfI    | 32 bp | fixed-byte |
| 16 | BslFI   | 25 bp | fixed-byte |

Enzymes can be specified by name (`--site BcgI`) or numeric ID (`--site 5`).

---

## Output Directory Structure

```
results/
├── 01_extract/                    # Step 1: Tag extraction
│   ├── sample1.BcgI.iibsp         # Binary tag file
│   ├── sample1.BcgI.stat.tsv      # Statistics
│   └── .done
│
├── 02_db_qual/                    # Step 2: Qualitative database
│   ├── BcgI.enzyme.iibdb          # All genome tags
│   ├── BcgI.species.iibdb         # Species-unique tags
│   ├── abfh_classify_with_speciename.txt.gz
│   └── .done
│
├── 02_db_quan/                    # Per-sample quantitative databases
│   ├── sample1/
│   │   ├── BcgI.species.iibdb
│   │   └── abfh_classify_with_speciename.txt.gz
│   └── sample2/
│
├── qualitative/                   # Qualitative screening results
│   ├── sample1/
│   │   ├── sample1.BcgI.xls
│   │   └── sample1.BcgI.GCF_detected.xls
│   └── .done
│
├── quantitative_sdb/              # Per-sample genome selection lists
│   ├── sample1/sdb.list
│   ├── sample2/sdb.list
│   └── .done
│
├── 04_quantify/                   # Quantitative profiling results
│   ├── sample1/
│   │   ├── sample1/sample1.BcgI.xls
│   │   └── .done
│   └── sample2/
│
└── 05_merge/                      # Final results
    ├── run1.all.xls               # Merged species abundance (all samples)
    ├── run1.filtered.xls          # Filtered (mock/control removed)
    ├── run1.func.xls              # Functional abundance (if --ko-mapping used)
    └── .done
```

---

## Binary File Format

Fast2bRAD-M uses a compact binary format (`.iibsp` / `.iibdb`) for storing hashed 2bRAD tags:

Two record streams share the same tooling (`inspect` auto-detects both):

- **Reference genome / database streams** (`.iibdb`): one record per tag,
  `[8-byte u64 hash][2-byte u16 id_length][id_bytes...]`. The id carries the
  contig (or `contig|offset` with `--record-pos`, or `gcf|idx|scaffold|pos|..`
  in a built database).
- **Sample tag streams** (`.iibsp`): an 8-byte header (`IIBS` + `u32` version)
  followed by bare `[8-byte u64 hash]` records. Read names are not stored —
  nothing downstream consumes them, and dropping them removes ~15 bytes per tag
  (hundreds of MB on a deeply sequenced sample). Older `.iibsp` files that do
  carry ids are still read correctly.

Tags are stored as canonical (lexicographically smaller of forward/reverse-complement)
FxHash values, which makes the format compact and fast to stream.

---

## Citation

If you use Fast2bRAD-M in your research, please cite the original 2bRAD-M paper:

> **2bRAD-M: Genome-level microbiome analysis using 2bRAD sequencing**

---

## License

Inherits the license of the original [2bRAD-M](https://github.com/HuangShiLab/2bRAD-M) project.
