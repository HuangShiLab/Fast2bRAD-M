# Fast2bRAD-M

**Fast2bRAD-M** is a high-performance Rust reimplementation of the [2bRAD-M](https://github.com/HuangShiLab/2bRAD-M) microbiome profiling pipeline. It delivers the same analytical results as the original Perl/Shell pipeline while achieving dramatically higher throughput through native parallelism and optimized I/O.

---

## Table of Contents

- [Features](#features)
- [Installation](#installation)
- [Quick Start](#quick-start)
- [Pipeline Overview](#pipeline-overview)
- [Subcommands](#subcommands)
  - [extract](#extract)
  - [build-qual-db](#build-qual-db)
  - [build-quan-db](#build-quan-db)
  - [quantify](#quantify)
  - [find-genome](#find-genome)
  - [merge](#merge)
  - [predict](#predict)
  - [classify](#classify)
  - [pipeline](#pipeline)
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
  -n 0.08 \                        # max N ratio
  -q 30 \                          # min quality score
  -p 80                            # min quality percent
```

**Input types**:

| Type | Description |
|------|-------------|
| 1    | Reference genome FASTA (for database construction) |
| 2    | Shotgun metagenome reads (SE or PE; PE can use PEAR merging) |
| 3    | Single 2bRAD tag reads |

**Output**:
- `{prefix}.{enzyme}.iibsp` — Binary tag file for sample reads (Types 2 & 3)
- `{prefix}.{enzyme}.iibdb` — Binary tag file for reference genomes (Type 1)
- `{prefix}.{enzyme}.stat.tsv` — Digest statistics

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

| ID | Name    | Tag Length |
|----|---------|-----------|
| 1  | CspCI   | 36 bp |
| 2  | AloI    | 37 bp |
| 3  | BsaXI   | 32 bp |
| 4  | BaeI    | 36 bp |
| **5**  | **BcgI** *(recommended)* | **32 bp** |
| 6  | CjeI    | 37 bp |
| 7  | PpiI    | 35 bp |
| 8  | PsrI    | 35 bp |
| 9  | BplI    | 35 bp |
| 10 | FalI    | 36 bp |
| 11 | Bsp24I  | 36 bp |
| 12 | HaeIV   | 37 bp |
| 13 | CjePI   | 38 bp |
| 14 | Hin4I   | 35 bp |
| 15 | AlfI    | 33 bp |
| 16 | BslFI   | 33 bp |

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

- Each record: `[8-byte u64 hash][4-byte u32 id_length][id_bytes...]`
- Tags are stored as canonical (lexicographically smaller of forward/reverse-complement) FxHash values
- This format enables fast random-access loading and minimal I/O

---

## Citation

If you use Fast2bRAD-M in your research, please cite the original 2bRAD-M paper:

> **2bRAD-M: Genome-level microbiome analysis using 2bRAD sequencing**

---

## License

Inherits the license of the original [2bRAD-M](https://github.com/HuangShiLab/2bRAD-M) project.
