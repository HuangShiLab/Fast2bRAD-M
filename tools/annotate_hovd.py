#!/usr/bin/env python3
"""
annotate_hovd.py

对 quantify_vip2b.py 使用 HOVD 数据库产生的合并病毒丰度表进行注释。
HOVD 的 metadata 字段与 VIP2B 不同，因此使用独立脚本。

输入:
  -i  合并丰度表 (fast2bRAD-M merge 输出的 *.all.xls)
  -d  HOVD metadata.tsv.gz（由 build_hovd_db.py 生成）
  -o  输出目录

输出:
  {outdir}/metadata_stats.tsv  — 按样品汇总的 contig metadata
  {outdir}/viral_taxonomy/{kingdom..species}_abund.tsv

用法:
  python tools/annotate_hovd.py \
    -i hovd_results/HOVD.all.xls \
    -d hovd_db/metadata.tsv.gz \
    -o hovd_results/annotation/
"""
import argparse
import gzip
import sys
from pathlib import Path
from typing import Dict, List, Tuple


def load_hovd_metadata(path: Path) -> Dict[str, Dict[str, str]]:
    """返回 {contig_id: {col: value}}。"""
    opener = gzip.open if str(path).endswith(".gz") else open
    meta: Dict[str, Dict[str, str]] = {}
    with opener(path, "rt", encoding="utf-8") as fh:
        header = fh.readline().rstrip("\n\r").split("\t")
        for line in fh:
            line = line.rstrip("\n\r")
            if not line:
                continue
            parts = line.split("\t")
            if len(parts) < len(header):
                continue
            cid = parts[0].strip()
            if not cid:
                continue
            meta[cid] = {h: parts[i].strip() for i, h in enumerate(header)}
    return meta


def parse_abundance(path: Path) -> Tuple[List[str], List[str], List[Dict[str, object]]]:
    """解析 merge 输出的 .all.xls。

    返回 (taxonomy_cols, sample_cols, rows)。
    """
    opener = gzip.open if str(path).endswith(".gz") else open
    rows: List[Dict[str, object]] = []
    with opener(path, "rt", encoding="utf-8") as fh:
        header = fh.readline().rstrip("\n\r").split("\t")
        taxonomy_cols = []
        sample_cols = []
        for h in header:
            if h in ("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species", "Strain"):
                taxonomy_cols.append(h)
            else:
                sample_cols.append(h)
        species_idx = taxonomy_cols.index("Species")
        for line in fh:
            line = line.rstrip("\n\r")
            if not line:
                continue
            parts = line.split("\t")
            if len(parts) < len(header):
                continue
            taxonomy = parts[: len(taxonomy_cols)]
            cid = taxonomy[species_idx].strip()
            abund = {}
            for j, sample in enumerate(sample_cols):
                val = parts[len(taxonomy_cols) + j]
                try:
                    abund[sample] = float(val)
                except ValueError:
                    abund[sample] = 0.0
            rows.append({"taxonomy": taxonomy, "cid": cid, "abund": abund})
    return taxonomy_cols, sample_cols, rows


def generate_metadata_stats(
    rows: List[Dict[str, object]],
    meta: Dict[str, Dict[str, str]],
    sample_cols: List[str],
    outdir: Path,
):
    """汇总 metadata 中的分类/质量/地理/宿主来源信息。"""
    # 选取有意义且数值/类别可汇总的字段
    categorical_fields = ["type", "checkv_quality", "provirus", "country", "continents", "site"]

    result: Dict[str, Dict[str, float]] = {}
    for field in categorical_fields:
        result[field] = {}

    for row in rows:
        cid = row["cid"]
        record = meta.get(cid, {})
        abund = row["abund"]
        for field in categorical_fields:
            val = record.get(field, "unknown").strip() or "unknown"
            key = f"{field}:{val}"
            d = result[field].setdefault(key, {s: 0.0 for s in sample_cols})
            for s in sample_cols:
                d[s] += abund[s]

    outpath = outdir / "metadata_stats.tsv"
    with open(outpath, "w") as fh:
        fh.write("Metric\t" + "\t".join(sample_cols) + "\n")
        for field in categorical_fields:
            for key in sorted(result[field].keys()):
                vals = [f"{result[field][key][s]:.6f}" for s in sample_cols]
                fh.write(f"{key}\t" + "\t".join(vals) + "\n")
    sys.stderr.write(f"Wrote {outpath}\n")


def generate_viral_taxonomy(
    rows: List[Dict[str, object]],
    sample_cols: List[str],
    outdir: Path,
):
    """按 Kingdom..Species 7 级汇总丰度。"""
    levels = ["kingdom", "phylum", "class", "order", "family", "genus", "species"]
    tax_dir = outdir / "viral_taxonomy"
    tax_dir.mkdir(parents=True, exist_ok=True)

    for idx, level_name in enumerate(levels):
        sums: Dict[str, Dict[str, float]] = {}
        for row in rows:
            tax = row["taxonomy"][: idx + 1]
            path = ";".join(tax) if all(tax) else "unknown"
            d = sums.setdefault(path, {s: 0.0 for s in sample_cols})
            for s in sample_cols:
                d[s] += row["abund"][s]

        outpath = tax_dir / f"{level_name}_abund.tsv"
        with open(outpath, "w") as fh:
            fh.write("Taxonomy\t" + "\t".join(sample_cols) + "\n")
            for path in sorted(sums.keys()):
                vals = [f"{sums[path][s]:.6f}" for s in sample_cols]
                fh.write(f"{path}\t" + "\t".join(vals) + "\n")
        sys.stderr.write(f"Wrote {outpath}\n")


def main():
    ap = argparse.ArgumentParser(description="Annotate HOVD viral abundance table with metadata")
    ap.add_argument("-i", "--input", required=True, help="Merged abundance table (*.all.xls)")
    ap.add_argument("-d", "--metadata", required=True, help="HOVD metadata.tsv.gz")
    ap.add_argument("-o", "--outdir", required=True, help="Output directory")
    args = ap.parse_args()

    inpath = Path(args.input)
    metapath = Path(args.metadata)
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    if not inpath.exists():
        sys.stderr.write(f"ERROR: input not found: {inpath}\n")
        sys.exit(1)
    if not metapath.exists():
        sys.stderr.write(f"ERROR: metadata not found: {metapath}\n")
        sys.exit(1)

    sys.stderr.write("Loading HOVD metadata ...\n")
    meta = load_hovd_metadata(metapath)
    sys.stderr.write(f"  {len(meta)} contigs loaded\n")

    sys.stderr.write("Loading abundance table ...\n")
    taxonomy_cols, sample_cols, rows = parse_abundance(inpath)
    sys.stderr.write(f"  {len(rows)} taxa, {len(sample_cols)} samples\n")

    generate_metadata_stats(rows, meta, sample_cols, outdir)
    generate_viral_taxonomy(rows, sample_cols, outdir)

    sys.stderr.write("All annotation tables done.\n")


if __name__ == "__main__":
    main()
