# Fast2bRAD-M viral reference databases

Fast2bRAD-M supports reference-based viral profiling with eight Type IIB
enzymes:

`AlfI`, `BcgI`, `BslFI`, `CjeI`, `CjePI`, `FalI`, `HaeIV`, `Hin4I`

## Naming: VIP2B vs UHGV

**VIP2B** (VIrome Profiler with type IIB restriction sites) is profiling
software. It is not the name of the viral reference catalog.

The gut-virome reference used by VIP2B is **UHGV** (Unified Human Gut Virome).
For Fast2bRAD-M, the converted distributable index is called **UHGV-8E**:

> **UHGV-8E** = the UHGV-based, eight-enzyme species index distributed with
> VIP2B v1.1 and converted to Fast2bRAD-M `CompactDatabase` format.

The historical wrapper names `convert_vip2b_db.py`, `quantify_vip2b.py`, and
`annotate_vip2b.py` remain for compatibility. New commands and output prefixes
should use `UHGV`.

## Available databases

| Database | Source | Scope | Entries | Compact tag records |
|---|---|---|---:|---:|
| **UHGV-8E** | [UHGV](https://uhgv.jgi.doe.gov/) via [VIP2B v1.1](https://zenodo.org/records/18944630) | Human gut virome / phage | 204,364 internal genome entries; 57,231 unique vOTU-level labels | 22,043,368 |
| **HOVD** | [HOVD](https://hovd.org/) | Human oral virome (OPD + OED) | 24,523 contigs | 7,747,787 |

## UHGV-8E

### Build from the VIP2B v1.1 source package

Download the three files from the
[VIP2B v1.1 Zenodo record](https://zenodo.org/records/18944630):

| Source file | MD5 |
|---|---|
| `8Enzyme.Species.uniq.marisa` | `0befbe79e1de76049f43babd33478dad` |
| `abfh_classify_with_speciename.txt.gz` | `282a34d153bee2a1f0000423cbee9d64` |
| `metadata.tsv.gz` | `f539c58de445b388bdf0457a412b77d3` |

Convert them to Fast2bRAD-M format:

```bash
python tools/convert_vip2b_db.py \
  -m 8Enzyme.Species.uniq.marisa \
  -c abfh_classify_with_speciename.txt.gz \
  -o uhgv_db/ \
  -l species \
  -s BcgI
```

`BcgI` is only a placeholder in the output filename. The resulting database
contains tags from all eight enzymes listed above.

### Converted package contents

| File | Size | MD5 |
|---|---:|---|
| `BcgI.species.iibdb` | 150,709,946 bytes | `e90a1295c7af10f1f39ab66138b26750` |
| `abfh_classify_with_speciename.txt.gz` | 987,119 bytes | `36ea92e27870025f98d5ab471b7a4e1f` |
| `metadata.tsv.gz` | 30,161,540 bytes | `f539c58de445b388bdf0457a412b77d3` |

Build statistics:

| Metric | Value |
|---|---:|
| Input marisa keys | 44,340,782 |
| Internal genome entries retained | 204,364 |
| Compact tag records | 22,043,368 |
| Unique vOTU-level labels in classify table | 57,231 |
| Metadata records | 168,536 |

The metadata includes `uhgv_genome`, `uhgv_votu`, lifestyle, host taxonomy,
viral taxonomy, UniRef90 annotations, and info-cluster annotations.

## HOVD

Download from [HOVD](https://hovd.org/):

- `HOVD-geneseqences.fasta`
- `HOVD-annotations.xlsx`

Build the Fast2bRAD-M database:

```bash
python tools/build_hovd_db.py \
  --fasta HOVD-geneseqences.fasta \
  --annotation HOVD-annotations.xlsx \
  -o hovd_db/ \
  -j 8
```

HOVD marks many OPD contigs as unclassified at genus/species levels. The builder
converts `uc_*` labels to `unknown` and uses the HOVD `contig_id` as the
species-level label so each abundance row remains identifiable.

## Profile samples

Prepare `samples.tsv`:

```tsv
sample1<TAB>/path/sample1_R1.fq.gz<TAB>/path/sample1_R2.fq.gz
sample2<TAB>/path/sample2_R1.fq.gz
```

Profile with UHGV-8E:

```bash
python tools/quantify_vip2b.py \
  -i samples.tsv \
  -d uhgv_db/ \
  -o uhgv_results/ \
  -j 16 \
  --qc no \
  --merge-prefix UHGV
```

Annotate UHGV-8E results:

```bash
python tools/annotate_vip2b.py \
  -i uhgv_results/UHGV.all.xls \
  -d uhgv_db/metadata.tsv.gz \
  -o uhgv_results/annotation/
```

Profile with HOVD:

```bash
python tools/quantify_vip2b.py \
  -i samples.tsv \
  -d hovd_db/ \
  -o hovd_results/ \
  -j 16 \
  --qc no \
  --merge-prefix HOVD

python tools/annotate_hovd.py \
  -i hovd_results/HOVD.all.xls \
  -d hovd_db/metadata.tsv.gz \
  -o hovd_results/annotation/
```

For already-demultiplexed 2bRAD tag reads, add `--input-type 3`. For WGS or
shotgun reads, use `--input-type 2`.

## Release artifact

The canonical Fast2bRAD-M converted package name is
`fast2bradm-uhgv-8e-v1.1`. A GitHub release asset should contain the four files
in the converted-package table plus a `CHECKSUMS.txt` file. Do not commit the
database binaries directly to git; the compact UHGV-8E index is approximately
144 MB.
