#!/usr/bin/env python3
"""
annotate_vip2b.py

对 quantify_vip2b.py 产生的合并病毒丰度表（如 VIP2B.all.xls）结合
metadata.tsv.gz 进行病毒注释：生活方式、宿主taxonomy、病毒taxonomy、
UniRef90 功能基因和 info 聚类。

输入:
  -i  合并丰度表 (fast2bRAD-M merge 输出的 *.all.xls)
  -d  VIP2B metadata.tsv.gz
  -o  输出目录

输出:
  {outdir}/Phenotype.tsv
  {outdir}/viral_taxonomy/{kingdom,phylum,class,order,family,genus}_abund.tsv
  {outdir}/host_taxonomy/{domain,phylum,class,order,family,genus,species}_abund.tsv
  {outdir}/viral_function/Uniref90.tsv
  {outdir}/viral_function/cluster.tsv

用法:
  python tools/annotate_vip2b.py \
    -i vip2b_results/VIP2B.all.xls \
    -d vip2b_db/metadata.tsv.gz \
    -o vip2b_results/annotation/
"""
import argparse
import gzip
import os
import sys
from pathlib import Path
from typing import Dict, List, Optional, Tuple


def load_metadata(path: Path) -> Dict[str, Dict[str, str]]:
    """返回 {vOTU: {col: value}}。metadata 第一列为 uhgv_genome，第二列为 uhgv_votu。"""
    opener = gzip.open if str(path).endswith(".gz") else open
    meta: Dict[str, Dict[str, str]] = {}
    with opener(path, "rt", encoding="utf-8") as fh:
        header = fh.readline().rstrip("\n\r").split("\t")
        if not header or header[0] != "uhgv_genome":
            sys.stderr.write("ERROR: metadata.tsv.gz header missing uhgv_genome\n")
            sys.exit(1)
        votu_idx = header.index("uhgv_votu") if "uhgv_votu" in header else 1
        for line in fh:
            line = line.rstrip("\n\r")
            if not line:
                continue
            parts = line.split("\t")
            if len(parts) < len(header):
                continue
            votu = parts[votu_idx].strip()
            if not votu:
                continue
            record = {h: parts[i].strip() for i, h in enumerate(header)}
            # 同 vOTU 可能对应多个 genome；保留第一次出现的信息
            if votu not in meta:
                meta[votu] = record
    return meta


def parse_abundance(path: Path) -> Tuple[List[str], List[str], List[Dict[str, object]]]:
    """解析 merge 输出的 .all.xls。

    返回 (taxonomy_cols, sample_cols, rows)，其中 rows 每项是
    {'taxonomy': [...], 'votu': str, 'abund': {sample: float}}。
    """
    opener = gzip.open if str(path).endswith(".gz") else open
    rows: List[Dict[str, object]] = []
    with opener(path, "rt", encoding="utf-8") as fh:
        header = fh.readline().rstrip("\n\r").split("\t")
        # fast2bRAD-M merge 的表头：Kingdom Phylum Class Order Family Genus Species sample1 ...
        taxonomy_cols = []
        sample_cols = []
        for i, h in enumerate(header):
            if h in ("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species", "Strain"):
                taxonomy_cols.append(h)
            else:
                sample_cols.append(h)
        if "Species" not in taxonomy_cols:
            sys.stderr.write(
                "ERROR: abundance table must contain Species column (VIP2B vOTU)\n"
            )
            sys.exit(1)
        species_idx = taxonomy_cols.index("Species")
        for line in fh:
            line = line.rstrip("\n\r")
            if not line:
                continue
            parts = line.split("\t")
            if len(parts) < len(header):
                continue
            taxonomy = parts[: len(taxonomy_cols)]
            votu = taxonomy[species_idx].strip()
            abund = {}
            for j, sample in enumerate(sample_cols):
                val = parts[len(taxonomy_cols) + j]
                try:
                    abund[sample] = float(val)
                except ValueError:
                    abund[sample] = 0.0
            rows.append({"taxonomy": taxonomy, "votu": votu, "abund": abund})
    return taxonomy_cols, sample_cols, rows


def safe_split(s: Optional[str], sep: str = "___") -> List[str]:
    if not s or s in ("-", "unknown"):
        return []
    return [x.strip() for x in s.split(sep) if x.strip()]


def generate_phenotype(
    rows: List[Dict[str, object]],
    meta: Dict[str, Dict[str, str]],
    sample_cols: List[str],
    outdir: Path,
):
    metrics = {
        "is_jumbo_phage": ["Yes", "No", "unknown"],
        "lifestyle": ["lytic", "temperate", "unknown"],
        "viralverify_prediction": ["Chromosome", "Plasmid", "Uncertain", "Virus", "unknown"],
    }
    result: Dict[str, Dict[str, float]] = {}
    for metric, values in metrics.items():
        for val in values:
            result[f"{metric}:{val}"] = {s: 0.0 for s in sample_cols}

    for row in rows:
        votu = row["votu"]
        record = meta.get(votu, {})
        abund = row["abund"]
        for metric, values in metrics.items():
            raw = record.get(metric, "unknown").strip()
            if raw not in values:
                raw = "unknown"
            key = f"{metric}:{raw}"
            for s in sample_cols:
                result[key][s] += abund[s]

    outpath = outdir / "Phenotype.tsv"
    with open(outpath, "w") as fh:
        fh.write("Metric\t" + "\t".join(sample_cols) + "\n")
        for key in sorted(result.keys()):
            vals = [f"{result[key][s]:.6f}" for s in sample_cols]
            fh.write(f"{key}\t" + "\t".join(vals) + "\n")
    sys.stderr.write(f"Wrote {outpath}\n")


def generate_functional(
    rows: List[Dict[str, object]],
    meta: Dict[str, Dict[str, str]],
    sample_cols: List[str],
    outdir: Path,
):
    uniref_sum: Dict[str, Dict[str, float]] = {}
    cluster_sum: Dict[str, Dict[str, float]] = {}

    for row in rows:
        votu = row["votu"]
        record = meta.get(votu, {})
        abund = row["abund"]
        for gene in safe_split(record.get("UniRef90_annotations")):
            d = uniref_sum.setdefault(gene, {s: 0.0 for s in sample_cols})
            for s in sample_cols:
                d[s] += abund[s]
        for clu in safe_split(record.get("info_annotations")):
            d = cluster_sum.setdefault(clu, {s: 0.0 for s in sample_cols})
            for s in sample_cols:
                d[s] += abund[s]

    func_dir = outdir / "viral_function"
    func_dir.mkdir(parents=True, exist_ok=True)

    def write_table(data: Dict[str, Dict[str, float]], outpath: Path, id_col: str):
        # 按每个 sample 的总和归一化
        totals = {s: sum(v[s] for v in data.values()) for s in sample_cols}
        with open(outpath, "w") as fh:
            fh.write(f"{id_col}\t" + "\t".join(sample_cols) + "\n")
            for key in sorted(data.keys()):
                vals = []
                for s in sample_cols:
                    total = totals[s]
                    v = data[key][s] / total if total > 0 else 0.0
                    vals.append(f"{v:.6f}")
                fh.write(f"{key}\t" + "\t".join(vals) + "\n")
        sys.stderr.write(f"Wrote {outpath}\n")

    write_table(uniref_sum, func_dir / "Uniref90.tsv", "gene")
    write_table(cluster_sum, func_dir / "cluster.tsv", "cluster")


def get_tax_level(tax_str: str, prefixes: List[str], target_idx: int) -> str:
    """从分号分隔的 taxonomy 字符串中提取到指定层级（含前缀）的路径。"""
    parts = [p.strip() for p in tax_str.split(";") if p.strip()]
    out = []
    for i in range(target_idx + 1):
        if i < len(parts):
            out.append(parts[i])
        else:
            out.append(f"{prefixes[i]}unknown")
    return ";".join(out)


def generate_taxonomy_breakdown(
    rows: List[Dict[str, object]],
    meta: Dict[str, Dict[str, str]],
    sample_cols: List[str],
    outdir: Path,
):
    configs = {
        "viral_taxonomy": {
            "prefixes": ["k__", "p__", "c__", "o__", "f__", "g__"],
            "levels": ["kingdom", "phylum", "class", "order", "family", "genus"],
            "save_dir": outdir / "viral_taxonomy",
        },
        "host_taxonomy": {
            "prefixes": ["d__", "p__", "c__", "o__", "f__", "g__", "s__"],
            "levels": ["domain", "phylum", "class", "order", "family", "genus", "species"],
            "save_dir": outdir / "host_taxonomy",
        },
    }

    for col, cfg in configs.items():
        cfg["save_dir"].mkdir(parents=True, exist_ok=True)
        for idx, level_name in enumerate(cfg["levels"]):
            sums: Dict[str, Dict[str, float]] = {}
            for row in rows:
                votu = row["votu"]
                tax_str = meta.get(votu, {}).get(col, "")
                if not tax_str:
                    continue
                path = get_tax_level(tax_str, cfg["prefixes"], idx)
                d = sums.setdefault(path, {s: 0.0 for s in sample_cols})
                for s in sample_cols:
                    d[s] += row["abund"][s]
            outpath = cfg["save_dir"] / f"{level_name}_abund.tsv"
            with open(outpath, "w") as fh:
                fh.write("Taxonomy\t" + "\t".join(sample_cols) + "\n")
                for path in sorted(sums.keys()):
                    vals = [f"{sums[path][s]:.6f}" for s in sample_cols]
                    fh.write(f"{path}\t" + "\t".join(vals) + "\n")
            sys.stderr.write(f"Wrote {outpath}\n")


def main():
    ap = argparse.ArgumentParser(
        description="Annotate VIP2B/Fast2bRAD-M viral abundance table with metadata"
    )
    ap.add_argument("-i", "--input", required=True, help="Merged abundance table (*.all.xls)")
    ap.add_argument("-d", "--metadata", required=True, help="VIP2B metadata.tsv.gz")
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

    sys.stderr.write("Loading metadata ...\n")
    meta = load_metadata(metapath)
    sys.stderr.write(f"  {len(meta)} vOTUs loaded\n")

    sys.stderr.write("Loading abundance table ...\n")
    taxonomy_cols, sample_cols, rows = parse_abundance(inpath)
    sys.stderr.write(f"  {len(rows)} taxa, {len(sample_cols)} samples\n")

    generate_phenotype(rows, meta, sample_cols, outdir)
    generate_functional(rows, meta, sample_cols, outdir)
    generate_taxonomy_breakdown(rows, meta, sample_cols, outdir)

    sys.stderr.write("All annotation tables done.\n")


if __name__ == "__main__":
    main()
