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
    - [quantify](#quantify)
    - [find-genome](#find-genome)
    - [merge](#merge)
    - [predict](#predict)
    - [classify](#classify)
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
- **ML Contamination Classification** — ONNX-based classification to detect contaminated taxa
- **Host Genotyping** — `f2brad-host` builds a host tag database and calls genotypes from 2bRAD reads
- **Holo-2bRAD Integration** — `f2brad-holo` performs one-pass joint host genotyping + microbial profiling with microbial cross-assignment masking
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
      │
      ▼ (optional, requires ONNX model)
[8] classify         →  per-sample classification with Prediction labels
```

---

## Subcommands

The project provides three command-line binaries:

| Binary | Purpose |
|--------|---------|
| `fast2bRAD-M` | Core microbiome profiling pipeline (extract → build-db → quantify → merge → predict/classify) |
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
  -d database_dir/ \     # directory with BcgI.species.iibdb + classify file
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

### `classify`

ML-based contamination classification using an ONNX model. Adds a `Prediction` column to the quantify output for each taxonomic entry.

**Features used** (4-dim input):
1. `ln(Sequenced_Tag_Num / Theoretical_Tag_Num)` — coverage ratio
2. `ln(G_score)` — combined breadth × depth
3. `ln(Sequenced_Reads_Num / Sequenced_Tag_Num)` — average depth
4. `ln(Theoretical_Reads / Total_Reads)` — theoretical abundance

```bash
fast2bRAD-M classify \
  -i 04_quantify/sample1/sample1.BcgI.xls \
  -m contamination_model.onnx \
  -o sample1.BcgI.classified.xls
```

**Parameters**:
| Parameter | Required | Description |
|-----------|----------|-------------|
| `-i` / `--input` | Yes | Input abundance table from `quantify` step |
| `-m` / `--model` | Yes | ONNX model file path |
| `-o` / `--output` | Yes | Output file path |

**Output**:
- Same TSV format as input with an additional `Prediction` column (integer label from the ONNX model)

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

**Output**:
- `genotypes.vcf` — host genotype calls
- `species_counts.tsv` — microbial taxon counts (when `--microbe-db-dir` is provided)
- `holo_classify.tsv` — read-classification summary (host/microbe/ambiguous fractions)
- `sample.iibsp.gz` — optional sample tag stream for downstream `fast2bRAD-M quantify` (with `--output-iibsp`)

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

├── classify/                      # ML classification results (optional)
│   ├── sample1.BcgI.classified.xls   # Per-sample with Prediction column
│   └── sample2.BcgI.classified.xls
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
