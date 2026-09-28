#!/usr/bin/env python3
"""
convert_vip2b_db.py

将 VIP2B 的 marisa 病毒数据库转换为 fast2bRAD-M 可读取的 CompactDatabase v3
(.iibdb) 及配套的 abfh_classify_with_speciename.txt.gz。

VIP2B 数据库结构（8Enzyme.Species.uniq.marisa）:
  - key: 酶切 tag，固定 40 字符，前导 '0' 为填充（zfill(40)）
  - value: 8 位 genome ID，如 b'00206754'

转换逻辑:
  1. 去掉 key 的前导 '0'，得到真实 tag 序列。
  2. 只保留 canonical（正/反向互补中字典序较小者）方向的 tag，避免
     VIP2B 同时存入正/反向互补导致理论 tag 数被重复计算。
  3. 用 FxHasher64(canonical) 计算 hash，写入 .iibdb。
  4. 从 VIP2B 的 abfh_classify_with_speciename.txt.gz 提取物种注释，
     输出同名分类映射文件。

输出:
  {outdir}/{site}.{level}.iibdb
  {outdir}/abfh_classify_with_speciename.txt.gz
  {outdir}/{site}.{level}.iibdb.stats.txt

用法:
  python tools/convert_vip2b_db.py \
    -m 8Enzyme.Species.uniq.marisa \
    -c abfh_classify_with_speciename.txt.gz \
    -o vip2b_db/ \
    -l species \
    -s BcgI

注意:
  VIP2B 是 8 种酶（AlfI/BcgI/BslFI/CjeI/CjePI/FalI/HaeIV/Hin4I）合并的
  数据库。site 参数仅用于命名输出文件并作为 fast2bRAD-M quantify 的占位
  酶名，必须是一个 fast2bRAD-M 已识别的酶（默认 BcgI）。
"""
import argparse
import gzip
import struct
import sys
from pathlib import Path
from typing import Dict, List, Set

K = 0x517CC1B727220A95  # rustc_hash FxHasher constant
COMPACT_MAGIC = b"IIBC"
COMPACT_VERSION = 3

# C-level reverse-complement translation table (bytes)
_RC_TABLE = bytes.maketrans(b"ACGTN", b"TGCAN")


class FxHasher64:
    """与 Rust rustc_hash::FxHasher (64-bit) 行为一致的 Python 实现。"""

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
            val = int.from_bytes(ba[:8], "little")
            self.add_to_hash(val)
            ba = ba[8:]
        if len(ba) >= 4:
            val = int.from_bytes(ba[:4], "little")
            self.add_to_hash(val)
            ba = ba[4:]
        for byte in ba:
            self.add_to_hash(byte)

    def finish(self) -> int:
        return self.hash


def canonical_hash(seq: bytes) -> int:
    """计算 canonical 方向的 FxHasher64 hash。seq 必须已为大写 ACGT。"""
    rc = seq.translate(_RC_TABLE)[::-1]
    canonical = seq if seq <= rc else rc
    h = FxHasher64()
    h.write(canonical)
    return h.finish()


def read_classify(path: Path) -> Dict[str, List[str]]:
    """读取 abfh_classify_with_speciename.txt.gz，返回 {genome_id: [ranks]}。"""
    classify: Dict[str, List[str]] = {}
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt", encoding="utf-8") as fh:
        for line in fh:
            line = line.rstrip("\n\r")
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t")
            if len(parts) < 8:
                continue
            gid = parts[0].strip()
            classify[gid] = [p.strip() for p in parts[1:8]]
    return classify


def load_marisa(marisa_path: Path):
    try:
        import marisa_trie
    except ImportError as e:
        sys.stderr.write(
            "ERROR: marisa_trie is required. Install with: pip install marisa-trie\n"
        )
        raise SystemExit(1) from e

    t = marisa_trie.BytesTrie()
    t.load(str(marisa_path))
    return t


def iter_genome_ids(trie) -> Set[str]:
    """快速扫描 marisa value，返回所有出现的 genome ID。"""
    observed: Set[str] = set()
    for _, value in trie.items():
        # value 可能为单个 genome ID 或多个 ID 的拼接（每个 8 字节）
        for i in range(0, len(value), 8):
            observed.add(value[i : i + 8].decode("ascii", errors="replace"))
    return observed


def convert_vip2b_db(
    marisa_path: Path,
    classify_path: Path,
    outdir: Path,
    level: str,
    site: str,
):
    import zstandard as zstd

    outdir.mkdir(parents=True, exist_ok=True)
    classify = read_classify(classify_path)
    sys.stderr.write(f"Classify entries: {len(classify)}\n")

    trie = load_marisa(marisa_path)
    sys.stderr.write(f"Loaded marisa trie: {len(trie)} keys\n")

    # Pass 1: 仅收集 genome IDs，不做 canonical 判断
    sys.stderr.write("Pass 1: collecting genome IDs ...\n")
    observed_gids = iter_genome_ids(trie)
    # 只保留有分类注释的 genome
    observed_gids = {gid for gid in observed_gids if gid in classify}
    sys.stderr.write(f"  observed genomes with taxonomy: {len(observed_gids)}\n")

    gcf_order: List[str] = sorted(observed_gids)
    gcf_to_idx: Dict[str, int] = {gid: i for i, gid in enumerate(gcf_order)}

    # 输出分类映射文件
    classify_out = outdir / "abfh_classify_with_speciename.txt.gz"
    with gzip.open(classify_out, "wt", encoding="utf-8") as fh:
        for gid in gcf_order:
            fh.write("\t".join([gid] + classify[gid]) + "\n")
    sys.stderr.write(f"Wrote classify file: {classify_out}\n")

    # Pass 2: 写入 iibdb
    db_path = outdir / f"{site}.{level}.iibdb"
    records_written = 0
    skipped_no_idx = 0
    skipped_non_canonical = 0

    with open(db_path, "wb") as fh:
        fh.write(COMPACT_MAGIC)
        fh.write(struct.pack("<I", COMPACT_VERSION))
        fh.write(struct.pack("<Q", 0))  # record_count placeholder
        fh.write(struct.pack("<I", len(gcf_order)))
        for gid in gcf_order:
            b = gid.encode("utf-8")
            if len(b) > 65535:
                b = b[:65535]
            fh.write(struct.pack("<H", len(b)))
            fh.write(b)

        compressor = zstd.ZstdCompressor(level=3)
        sys.stderr.write("Pass 2: writing compact database ...\n")
        with compressor.stream_writer(fh, closefd=False) as writer:
            for tag_str, value in trie.items():
                seq = tag_str.lstrip("0").encode("ascii")
                if not seq:
                    continue
                # 只保留 canonical 方向，避免正/反向互补重复
                if seq > seq.translate(_RC_TABLE)[::-1]:
                    skipped_non_canonical += 1
                    continue
                h = canonical_hash(seq)
                for i in range(0, len(value), 8):
                    gid = value[i : i + 8].decode("ascii", errors="replace")
                    idx = gcf_to_idx.get(gid)
                    if idx is None:
                        skipped_no_idx += 1
                        continue
                    writer.write(struct.pack("<Q", h))
                    writer.write(struct.pack("<I", idx))
                    records_written += 1
                if records_written % 5_000_000 == 0:
                    sys.stderr.write(f"  {records_written:,} records written\n")
                    sys.stderr.flush()

        fh.seek(8)
        fh.write(struct.pack("<Q", records_written))

    stats_path = Path(str(db_path) + ".stats.txt")
    with open(stats_path, "w") as fh:
        fh.write(f"input_keys\t{len(trie)}\n")
        fh.write(f"genomes\t{len(gcf_order)}\n")
        fh.write(f"records_written\t{records_written}\n")
        fh.write(f"skipped_non_canonical\t{skipped_non_canonical}\n")
        fh.write(f"skipped_no_idx\t{skipped_no_idx}\n")
        fh.write(f"level\t{level}\n")
        fh.write(f"site\t{site}\n")

    sys.stderr.write(
        f"Wrote {db_path}: {records_written:,} records, "
        f"{len(gcf_order)} genomes\n"
    )


def main():
    ap = argparse.ArgumentParser(
        description="Convert VIP2B marisa viral DB to fast2bRAD-M .iibdb"
    )
    ap.add_argument("-m", "--marisa", required=True, help="Input marisa DB file")
    ap.add_argument(
        "-c", "--classify", required=True, help="VIP2B abfh_classify_with_speciename.txt.gz"
    )
    ap.add_argument("-o", "--outdir", required=True, help="Output directory")
    ap.add_argument(
        "-l", "--level", default="species", help="Taxonomy level (species/genus/...)"
    )
    ap.add_argument(
        "-s",
        "--site",
        default="BcgI",
        help="Placeholder enzyme name for output files (must be a valid fast2bRAD-M enzyme)",
    )
    args = ap.parse_args()

    convert_vip2b_db(
        Path(args.marisa),
        Path(args.classify),
        Path(args.outdir),
        args.level,
        args.site,
    )


if __name__ == "__main__":
    main()
