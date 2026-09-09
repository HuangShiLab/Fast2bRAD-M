#!/usr/bin/env python3
"""
build_hovd_db.py

把 HOVD (Human Oral Virome Database) 的 genome FASTA 和注释表格转换成
fast2bRAD-M 可用的 8 酶合并 CompactDatabase v3 (.iibdb) 及配套文件。

输入:
  --fasta        HOVD-geneseqences.fasta（含 OPD + OED）
  --annotation   HOVD-annotations.xlsx（含 OPD、OED 两个 sheet）
  -o             输出目录

输出:
  {outdir}/{site}.species.iibdb
  {outdir}/abfh_classify_with_speciename.txt.gz
  {outdir}/metadata.tsv.gz
  {outdir}/{site}.species.iibdb.stats.txt

用法:
  python tools/build_hovd_db.py \
    --fasta ~/Downloads/HOVD/HOVD-geneseqences.fasta \
    --annotation ~/Downloads/HOVD/HOVD-annotations.xlsx \
    -o hovd_db/ -j 8
"""
import argparse
import gzip
import re
import struct
import sys
from collections import defaultdict
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import zstandard as zstd

K = 0x517CC1B727220A95  # rustc_hash FxHasher constant
COMPACT_MAGIC = b"IIBC"
COMPACT_VERSION = 3

# VIP2B 8 酶正则（与 2bRAD-M / VIP2B 一致）
ENZYME_PATTERNS = {
    "AlfI": r"(?=([AGCT]{10}GCA[AGCT]{6}TGC[AGCT]{10}))",
    "BcgI": r"(?=([AGCT]{10}CGA[AGCT]{6}TGC[AGCT]{10}))",
    "BslFI": r"(?=([AGCT]{6}GGGAC[AGCT]{14}))",
    "CjeI": r"(?=([AGCT]{8}CCA[AGCT]{6}GT[AGCT]{9}))",
    "CjePI": r"(?=([AGCT]{7}CCA[AGCT]{7}TC[AGCT]{8}))",
    "FalI": r"(?=([AGCT]{8}AAG[AGCT]{5}CTT[AGCT]{8}))",
    "HaeIV": r"(?=([AGCT]{7}GA[CT][AGCT]{5}[AG]TC[AGCT]{9}))",
    "Hin4I": r"(?=([AGCT]{8}GA[CT][AGCT]{5}[GAC]TC[AGCT]{8}))",
}

_RC_TABLE = bytes.maketrans(b"ACGTN", b"TGCAN")


class FxHasher64:
    __slots__ = ("hash",)

    def __init__(self):
        self.hash = 0

    def add_to_hash(self, i: int):
        h = self.hash
        h = ((h << 5) | (h >> 59)) & 0xFFFFFFFFFFFFFFFF
        h ^= i & 0xFFFFFFFFFFFFFFFF
        h = (h * K) & 0xFFFFFFFFFFFFFFFF
        self.hash = h

    def write(self, data: bytes):
        ba = bytearray(data)
        while len(ba) >= 8:
            self.add_to_hash(int.from_bytes(ba[:8], "little"))
            ba = ba[8:]
        if len(ba) >= 4:
            self.add_to_hash(int.from_bytes(ba[:4], "little"))
            ba = ba[4:]
        for byte in ba:
            self.add_to_hash(byte)

    def finish(self) -> int:
        return self.hash


def canonical_hash(seq: bytes) -> int:
    rc = seq.translate(_RC_TABLE)[::-1]
    canonical = seq if seq <= rc else rc
    h = FxHasher64()
    h.write(canonical)
    return h.finish()


def parse_hovd_taxonomy(tax_str: str) -> List[str]:
    """把 HOVD taxonomy 字符串（species→kingdom）转成 Kingdom→Species 列表。

    HOVD 用 'uc_<rank>'（如 uc_Caudoviricetes）标记未分类层级；这里统一替换为
    'unknown'，使输出与 fast2bRAD-M 的 unknown 约定一致。
    """
    if not tax_str or not isinstance(tax_str, str):
        return ["unknown"] * 7
    parts = [p.strip() for p in tax_str.split(";") if p.strip()]
    if len(parts) == 1 and parts[0].lower() == "unclassified":
        return ["unknown"] * 7
    # HOVD 的 uc_ 前缀表示 unclassified
    parts = ["unknown" if p.lower().startswith("uc_") else p for p in parts]
    # 补齐到 7 级
    while len(parts) < 7:
        parts.insert(0, "unknown")
    if len(parts) > 7:
        parts = parts[:7]
    # 当前顺序是 species→kingdom，需要反转
    parts.reverse()
    return parts


def parse_annotations(anno_path: Path) -> Tuple[Dict[str, List[str]], Dict[str, Dict[str, str]]]:
    """解析 HOVD-annotations.xlsx，返回 (taxonomy_map, metadata_map)。"""
    try:
        import pandas as pd
    except ImportError as e:
        sys.stderr.write("ERROR: pandas is required. Install with: pip install pandas openpyxl\n")
        raise SystemExit(1) from e

    tax_map: Dict[str, List[str]] = {}
    meta_map: Dict[str, Dict[str, str]] = {}

    xl = pd.ExcelFile(str(anno_path))
    sheets = xl.sheet_names
    sys.stderr.write(f"Annotation sheets: {sheets}\n")

    # OPD sheet
    if "OPD" in sheets:
        df = xl.parse("OPD")
        for _, row in df.iterrows():
            cid = str(row.get("contig_id", "")).strip()
            if not cid:
                continue
            tax_map[cid] = parse_hovd_taxonomy(str(row.get("taxonomy", "")))
            meta_map[cid] = {
                "contig_id": cid,
                "type": "phage",
                "country": str(row.get("country", "")),
                "continents": str(row.get("continents", "")),
                "site": str(row.get("site", "")),
                "contig_length": str(row.get("contig_length", "")),
                "checkv_quality": str(row.get("checkv_quality", "")),
                "completeness": str(row.get("completeness", "")),
                "contamination": str(row.get("contamination", "")),
                "provirus": str(row.get("provirus", "")),
                "gene_count": str(row.get("gene_count", "")),
                "host_genes": str(row.get("host_genes", "")),
                "taxonomy": str(row.get("taxonomy", "")),
            }

    # OED sheet
    if "OED" in sheets:
        df = xl.parse("OED")
        for _, row in df.iterrows():
            cid = str(row.get("contig_id", "")).strip()
            if not cid:
                continue
            family = str(row.get("Family.level.taxonomy", "")).strip()
            tax = ["unknown", "unknown", "unknown", "unknown", family or "unknown", "unknown", "unknown"]
            tax_map[cid] = tax
            meta_map[cid] = {
                "contig_id": cid,
                "type": "eukaryotic_virus",
                "country": str(row.get("country", "")),
                "continents": str(row.get("continents", "")),
                "site": str(row.get("site", "")),
                "contig_length": str(row.get("contig_length", "")),
                "family": family,
                "completeness": str(row.get("completeness", "")),
                "contamination": str(row.get("contamination", "")),
                "taxonomy": family,
            }

    return tax_map, meta_map


def parse_fasta_headers(fasta_path: Path) -> List[str]:
    """只扫描 FASTA header，返回按出现顺序的 contig_id 列表。"""
    ids: List[str] = []
    opener = gzip.open if str(fasta_path).endswith(".gz") else open
    with opener(fasta_path, "rt", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if line.startswith(">"):
                cid = line[1:].split()[0].strip()
                if cid:
                    ids.append(cid)
    return ids


def iter_fasta_records(fasta_path: Path):
    """逐条产生 (contig_id, sequence_upper_bytes)。"""
    opener = gzip.open if str(fasta_path).endswith(".gz") else open
    with opener(fasta_path, "rt", encoding="utf-8", errors="replace") as fh:
        cid = None
        seq_parts: List[str] = []
        for line in fh:
            line = line.rstrip("\n\r")
            if line.startswith(">"):
                if cid is not None:
                    yield cid, "".join(seq_parts).upper().encode("ascii")
                cid = line[1:].split()[0].strip()
                seq_parts = []
            else:
                seq_parts.append(line)
        if cid is not None:
            yield cid, "".join(seq_parts).upper().encode("ascii")


def extract_tags(seq: bytes, pattern: re.Pattern) -> List[bytes]:
    """用 VIP2B 方式从正、反两条链提取 tag，返回 tag bytes 列表。"""
    seq_str = seq.decode("ascii")
    rc_str = seq.translate(_RC_TABLE)[::-1].decode("ascii")
    tags = []
    for s in (seq_str, rc_str):
        for m in pattern.finditer(s):
            tags.append(m.group(1).encode("ascii"))
    return tags


def build_hovd_db(
    fasta_path: Path,
    anno_path: Path,
    outdir: Path,
    site: str,
    enzymes: List[str],
    threads: int,
):
    outdir.mkdir(parents=True, exist_ok=True)
    tax_map, meta_map = parse_annotations(anno_path)
    sys.stderr.write(f"Annotation entries: {len(tax_map)}\n")

    # 只保留有注释的 contig
    fasta_ids = [cid for cid in parse_fasta_headers(fasta_path) if cid in tax_map]
    sys.stderr.write(f"FASTA contigs with annotation: {len(fasta_ids)}\n")
    if not fasta_ids:
        sys.stderr.write("ERROR: no matching contigs found\n")
        sys.exit(1)

    gcf_order = sorted(fasta_ids)
    gcf_to_idx = {cid: i for i, cid in enumerate(gcf_order)}

    # 写分类文件
    classify_out = outdir / "abfh_classify_with_speciename.txt.gz"
    with gzip.open(classify_out, "wt", encoding="utf-8") as fh:
        for cid in gcf_order:
            fh.write("\t".join([cid] + tax_map[cid]) + "\n")
    sys.stderr.write(f"Wrote classify file: {classify_out}\n")

    # 写 metadata 文件
    metadata_out = outdir / "metadata.tsv.gz"
    meta_cols = list(next(iter(meta_map.values())).keys())
    with gzip.open(metadata_out, "wt", encoding="utf-8") as fh:
        fh.write("\t".join(meta_cols) + "\n")
        for cid in gcf_order:
            fh.write("\t".join(str(meta_map[cid].get(c, "")) for c in meta_cols) + "\n")
    sys.stderr.write(f"Wrote metadata file: {metadata_out}\n")

    # 编译正则
    compiled = {e: re.compile(ENZYME_PATTERNS[e]) for e in enzymes}

    # 写 iibdb
    db_path = outdir / f"{site}.species.iibdb"
    records_written = 0
    skipped_no_anno = 0

    with open(db_path, "wb") as fh:
        fh.write(COMPACT_MAGIC)
        fh.write(struct.pack("<I", COMPACT_VERSION))
        fh.write(struct.pack("<Q", 0))
        fh.write(struct.pack("<I", len(gcf_order)))
        for cid in gcf_order:
            b = cid.encode("utf-8")
            if len(b) > 65535:
                b = b[:65535]
            fh.write(struct.pack("<H", len(b)))
            fh.write(b)

        compressor = zstd.ZstdCompressor(level=3)
        with compressor.stream_writer(fh, closefd=False) as writer:
            for cid, seq in iter_fasta_records(fasta_path):
                idx = gcf_to_idx.get(cid)
                if idx is None:
                    skipped_no_anno += 1
                    continue
                seen: set = set()
                for enzyme in enzymes:
                    for tag in extract_tags(seq, compiled[enzyme]):
                        h = canonical_hash(tag)
                        if h in seen:
                            continue
                        seen.add(h)
                        writer.write(struct.pack("<Q", h))
                        writer.write(struct.pack("<I", idx))
                        records_written += 1
                if records_written % 1_000_000 == 0:
                    sys.stderr.write(f"  {records_written:,} records written\n")
                    sys.stderr.flush()

        fh.seek(8)
        fh.write(struct.pack("<Q", records_written))

    stats_path = Path(str(db_path) + ".stats.txt")
    with open(stats_path, "w") as fh:
        fh.write(f"genomes\t{len(gcf_order)}\n")
        fh.write(f"records_written\t{records_written}\n")
        fh.write(f"skipped_no_anno\t{skipped_no_anno}\n")
        fh.write(f"enzymes\t{','.join(enzymes)}\n")
        fh.write(f"site\t{site}\n")

    sys.stderr.write(
        f"Wrote {db_path}: {records_written:,} records, {len(gcf_order)} genomes\n"
    )


def main():
    ap = argparse.ArgumentParser(description="Build HOVD 8-enzyme viral database for Fast2bRAD-M")
    ap.add_argument("--fasta", required=True, help="HOVD-geneseqences.fasta")
    ap.add_argument("--annotation", required=True, help="HOVD-annotations.xlsx")
    ap.add_argument("-o", "--outdir", required=True, help="Output directory")
    ap.add_argument("-s", "--site", default="BcgI", help="Placeholder enzyme name")
    ap.add_argument("-e", "--enzymes", default=",".join(ENZYME_PATTERNS.keys()), help="Comma-separated enzymes")
    ap.add_argument("-j", "--threads", type=int, default=4, help="Reserved for future parallelism")
    args = ap.parse_args()

    build_hovd_db(
        Path(args.fasta),
        Path(args.annotation),
        Path(args.outdir),
        args.site,
        [e.strip() for e in args.enzymes.split(",") if e.strip()],
        args.threads,
    )


if __name__ == "__main__":
    main()
