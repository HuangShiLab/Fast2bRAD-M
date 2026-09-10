#!/usr/bin/env python3
"""
compare_2brad_wms.py

比较同一样品（subject + timepoint）的 BcgI 2bRAD HOVD 病毒 profile 与
8 酶 WMS HOVD profile 的一致性。

输入:
  -2 2bRAD 合并表 (HOVD_2bRAD.all.xls)
  -w WMS 合并表 (HOVD_WMS.all.xls)
  -o 输出目录

输出:
  summary.tsv, per_pair_correlation.tsv, detection_overlap.tsv
"""
import argparse
import gzip
import sys
from pathlib import Path
from typing import Dict, List, Tuple


def parse_sample_name(name: str) -> Tuple[str, str]:
    """s11_12 -> (11, 12)"""
    name = name.lstrip("s")
    parts = name.rsplit("_", 1)
    return parts[0], parts[1]


def load_table(path: Path) -> Tuple[List[str], List[str], Dict[str, Dict[str, float]]]:
    """返回 (taxonomy_cols, sample_cols, {taxon_key: {sample: abund}})。"""
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt", encoding="utf-8") as fh:
        header = fh.readline().rstrip("\n\r").split("\t")
        taxonomy_cols = []
        sample_cols = []
        for h in header:
            hc = h.lstrip("#")
            if hc in ("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species", "Strain"):
                taxonomy_cols.append(hc)
            else:
                sample_cols.append(h)
        taxa_idx = {h: i for i, h in enumerate(header)}
        rows: Dict[str, Dict[str, float]] = {}
        for line in fh:
            parts = line.rstrip("\n\r").split("\t")
            if len(parts) < len(header):
                continue
            tax_key = "|".join(parts[tax_idx[h]] for h in header if h in taxa_idx and h.lstrip("#") in taxonomy_cols)
            # 用 Species 作为 taxon key 足够比较
            species_idx = header.index("Species") if "Species" in header else -1
            if species_idx >= 0:
                tax_key = parts[species_idx].strip()
            else:
                tax_key = "|".join(parts[:len(taxonomy_cols)])
            rows[tax_key] = {}
            for s in sample_cols:
                try:
                    rows[tax_key][s] = float(parts[tax_idx[s]])
                except (ValueError, KeyError):
                    rows[tax_key][s] = 0.0
        return taxonomy_cols, sample_cols, rows


def pearson_r(x: List[float], y: List[float]) -> float:
    n = len(x)
    if n == 0:
        return 0.0
    mx = sum(x) / n
    my = sum(y) / n
    num = sum((xi - mx) * (yi - my) for xi, yi in zip(x, y))
    denx = sum((xi - mx) ** 2 for xi in x) ** 0.5
    deny = sum((yi - my) ** 2 for yi in y) ** 0.5
    if denx == 0 or deny == 0:
        return 0.0
    return num / (denx * deny)


def bray_curtis(x: List[float], y: List[float]) -> float:
    s = sum(abs(xi - yi) for xi, yi in zip(x, y))
    t = sum(xi + yi for xi, yi in zip(x, y))
    if t == 0:
        return 1.0
    return s / t


def main():
    ap = argparse.ArgumentParser(description="Compare 2bRAD and WMS HOVD profiles")
    ap.add_argument("-2", "--rad", required=True, help="2bRAD merged abundance table")
    ap.add_argument("-w", "--wms", required=True, help="WMS merged abundance table")
    ap.add_argument("-o", "--outdir", required=True, help="Output directory")
    args = ap.parse_args()

    rad_path = Path(args.rad)
    wms_path = Path(args.wms)
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    _, rad_samples, rad_table = load_table(rad_path)
    _, wms_samples, wms_table = load_table(wms_path)

    # Match by subject_timepoint
    rad_lookup: Dict[Tuple[str, str], str] = {}
    for s in rad_samples:
        subj, tp = parse_sample_name(s)
        rad_lookup[(subj, tp)] = s

    pairs: List[Tuple[str, str, str, str]] = []  # (subject, timepoint, rad_sample, wms_sample)
    for ws in wms_samples:
        subj, tp = parse_sample_name(ws)
        if (subj, tp) in rad_lookup:
            pairs.append((subj, tp, rad_lookup[(subj, tp)], ws))

    if not pairs:
        sys.stderr.write("ERROR: no matched sample pairs found\n")
        sys.exit(1)

    sys.stderr.write(f"Matched {len(pairs)} sample pairs\n")

    # All common taxa
    all_taxa = sorted(set(rad_table.keys()) | set(wms_table.keys()))

    pair_rows = []
    for subj, tp, rs, ws in pairs:
        rad_vec = [rad_table.get(t, {}).get(rs, 0.0) for t in all_taxa]
        wms_vec = [wms_table.get(t, {}).get(ws, 0.0) for t in all_taxa]
        r = pearson_r(rad_vec, wms_vec)
        bc = bray_curtis(rad_vec, wms_vec)
        # detection overlap (taxa with >0 in both)
        rad_detected = {t for t in all_taxa if rad_table.get(t, {}).get(rs, 0.0) > 0}
        wms_detected = {t for t in all_taxa if wms_table.get(t, {}).get(ws, 0.0) > 0}
        overlap = len(rad_detected & wms_detected)
        jaccard = overlap / len(rad_detected | wms_detected) if (rad_detected | wms_detected) else 0.0
        pair_rows.append((subj, tp, rs, ws, r, 1 - bc, overlap, jaccard, len(rad_detected), len(wms_detected)))

    pair_path = outdir / "per_pair_correlation.tsv"
    with open(pair_path, "w") as fh:
        fh.write("subject\ttimepoint\trad_sample\twms_sample\tpearson_r\tbray_curtis_similarity\tdetection_overlap\tjaccard\trad_detected\twms_detected\n")
        for row in pair_rows:
            fh.write("\t".join(str(x) for x in row) + "\n")
    sys.stderr.write(f"Wrote {pair_path}\n")

    # Subject-level aggregation
    subject_groups: Dict[str, List[Tuple]] = {}
    for row in pair_rows:
        subject_groups.setdefault(row[0], []).append(row)

    summary_rows = []
    for subj, rows in sorted(subject_groups.items()):
        avg_r = sum(r[4] for r in rows) / len(rows)
        avg_bc = sum(r[5] for r in rows) / len(rows)
        avg_jaccard = sum(r[7] for r in rows) / len(rows)
        summary_rows.append((subj, len(rows), avg_r, avg_bc, avg_jaccard))

    summary_path = outdir / "summary.tsv"
    with open(summary_path, "w") as fh:
        fh.write("subject\tn_timepoints\tmean_pearson_r\tmean_bray_curtis_similarity\tmean_jaccard\n")
        for row in summary_rows:
            fh.write("\t".join(f"{x:.6f}" if isinstance(x, float) else str(x) for x in row) + "\n")
        # Overall
        all_r = [r[4] for r in pair_rows]
        all_bc = [r[5] for r in pair_rows]
        all_j = [r[7] for r in pair_rows]
        fh.write(f"overall\t{len(pair_rows)}\t{sum(all_r)/len(all_r):.6f}\t{sum(all_bc)/len(all_bc):.6f}\t{sum(all_j)/len(all_j):.6f}\n")
    sys.stderr.write(f"Wrote {summary_path}\n")


if __name__ == "__main__":
    main()
