#!/usr/bin/env python3
"""
convert_fasta_db.py

将 2bRAD-M 风格的预消化 tag FASTA (如 BcgI.species.fa.gz) 转换为
fast2bRAD-M f2brad-holo / f2brad-m quantify 使用的 CompactDatabase v3 (.iibdb)。

输入 FASTA 要求:
  - 每条记录是一个已消化的 tag (默认 32 bp)
  - header 第一个字段为 GCF/GCA ID, 如:
      >GCF_014075335.1|1|NZ_CP050047.1|353|1|0
  - 同文件内可包含多个基因组, 按 header 首字段分组

输出:
  - {out_prefix}.iibdb   CompactDatabase v3 (zstd 压缩 records)
  - {out_prefix}.iibdb.stats.txt  统计信息

用法:
  python tools/convert_fasta_db.py \
    -i /path/to/BcgI.species.fa.gz \
    -o /path/to/BcgI.species.iibdb \
    --enzyme BcgI
"""
import argparse
import gzip
import sys
import struct
from collections import defaultdict
from pathlib import Path
from typing import Dict, List, Tuple, Optional

K = 0x517cc1b727220a95  # rustc_hash FxHasher constant
COMPACT_MAGIC = b"IIBC"
COMPACT_VERSION = 3


class FxHasher64:
    """与 Rust rustc_hash::FxHasher (64-bit) 行为一致的 Python 实现。"""

    def __init__(self):
        self.hash = 0

    def add_to_hash(self, i: int):
        h = self.hash
        h = ((h << 5) | (h >> 59)) & 0xFFFFFFFFFFFFFFFF  # rotate_left(5)
        h ^= i & 0xFFFFFFFFFFFFFFFF
        h = (h * K) & 0xFFFFFFFFFFFFFFFF
        self.hash = h

    def write(self, data: bytes):
        """Matches fxhash 0.2.1 write64: 8-byte chunks, then a 4-byte chunk,
        then individual bytes."""
        ba = bytearray(data)
        # 64-bit native-endian chunks (x86_64 = little-endian)
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


def reverse_complement(seq: bytes) -> bytes:
    comp = bytearray(len(seq))
    for i, b in enumerate(seq):
        if b == 65 or b == 97:  # A/a
            comp[len(seq) - 1 - i] = 84  # T
        elif b == 84 or b == 116:  # T/t
            comp[len(seq) - 1 - i] = 65  # A
        elif b == 67 or b == 99:  # C/c
            comp[len(seq) - 1 - i] = 71  # G
        elif b == 71 or b == 103:  # G/g
            comp[len(seq) - 1 - i] = 67  # C
        elif b == 78 or b == 110:  # N/n
            comp[len(seq) - 1 - i] = 78  # N
        else:
            comp[len(seq) - 1 - i] = b
    return bytes(comp)


def canonical_hash(seq: bytes) -> int:
    """计算 canonical (取正/反向互补中字典序较小者) 的 FxHasher64 hash。"""
    seq = bytes(seq).upper()
    rc = reverse_complement(seq)
    canonical = seq if seq <= rc else rc
    h = FxHasher64()
    h.write(canonical)
    return h.finish()


def parse_fasta_records(path: Path):
    """逐条产生 (header_line_without_gt, sequence_bytes)。"""
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rb") as fh:
        header = None
        seq_parts: List[bytes] = []
        for raw in fh:
            line = raw.rstrip(b"\n\r")
            if not line:
                continue
            if line.startswith(b">"):
                if header is not None:
                    yield header, b"".join(seq_parts)
                header = line[1:].decode("ascii", errors="replace")
                seq_parts = []
            else:
                seq_parts.append(line)
        if header is not None:
            yield header, b"".join(seq_parts)


def extract_gcf_id(header: str) -> str:
    """从 FASTA header 提取 GCF/GCA ID (首字段)。"""
    first = header.split("|")[0].strip()
    # 兼容文件名形式, 如 GCF_014075335.1_genomic -> GCF_014075335.1
    name = first.split("/")[-1]
    if "_genomic" in name:
        name = name[: name.find("_genomic")]
    if name.startswith("GCF_") or name.startswith("GCA_"):
        parts = name.split("_")
        if len(parts) >= 2:
            return f"{parts[0]}_{parts[1]}"
    return name


def convert_fasta_to_iibdb(
    input_path: Path,
    output_path: Path,
    enzyme: str,
    expected_length: Optional[int] = None,
):
    sys.stderr.write(f"Parsing {input_path} ...\n")

    # First pass: collect GCF IDs in order of first appearance
    gcf_order: List[str] = []
    gcf_to_idx: Dict[str, int] = {}
    gcf_tag_counts: Dict[str, int] = defaultdict(int)

    # Also collect per-GCF tag sequences to compute hashes in memory.
    # For very large DBs this could be streamed twice; BcgI.species.fa.gz is ~14G,
    # so two passes are acceptable and simpler than external sort.
    per_gcf_tags: Dict[str, List[bytes]] = defaultdict(list)

    total_records = 0
    length_issues = 0

    for header, seq in parse_fasta_records(input_path):
        total_records += 1
        if expected_length and len(seq) != expected_length:
            length_issues += 1
        gcf_id = extract_gcf_id(header)
        if gcf_id not in gcf_to_idx:
            gcf_to_idx[gcf_id] = len(gcf_order)
            gcf_order.append(gcf_id)
        per_gcf_tags[gcf_id].append(seq)
        gcf_tag_counts[gcf_id] += 1

        if total_records % 5_000_000 == 0:
            sys.stderr.write(f"  {total_records:,} records parsed\n")

    sys.stderr.write(
        f"Done. {total_records:,} records, {len(gcf_order)} genomes, "
        f"length_issues={length_issues}\n"
    )

    # Compute canonical hashes and deduplicate globally:
    # If the same canonical hash maps to multiple GCFs, keep the first occurrence
    # (original 2bRAD-M pipeline should have already removed cross-species shared tags,
    #  so collisions should be rare).
    hash_to_gcf: Dict[int, int] = {}
    collisions = 0
    records_written = 0

    sys.stderr.write("Computing hashes and writing compact database ...\n")
    output_path.parent.mkdir(parents=True, exist_ok=True)

    import zstandard as zstd  # local import so hash-only imports work without it

    with open(output_path, "wb") as fh:
        # Header (uncompressed)
        fh.write(COMPACT_MAGIC)
        fh.write(struct.pack("<I", COMPACT_VERSION))
        fh.write(struct.pack("<Q", 0))  # record_count placeholder
        fh.write(struct.pack("<I", len(gcf_order)))
        for gcf_id in gcf_order:
            b = gcf_id.encode("utf-8")
            if len(b) > 65535:
                b = b[:65535]
            fh.write(struct.pack("<H", len(b)))
            fh.write(b)

        # Records (zstd compressed)
        compressor = zstd.ZstdCompressor(level=3)
        with compressor.stream_writer(fh) as writer:
            for gcf_id in gcf_order:
                idx = gcf_to_idx[gcf_id]
                for seq in per_gcf_tags[gcf_id]:
                    h = canonical_hash(seq)
                    existing = hash_to_gcf.get(h)
                    if existing is None:
                        hash_to_gcf[h] = idx
                        writer.write(struct.pack("<Q", h))
                        writer.write(struct.pack("<I", idx))
                        records_written += 1
                    elif existing != idx:
                        collisions += 1

        # Patch record_count
        fh.seek(8)
        fh.write(struct.pack("<Q", records_written))

    stats_path = Path(str(output_path) + ".stats.txt")
    with open(stats_path, "w") as fh:
        fh.write(f"input_records\t{total_records}\n")
        fh.write(f"genomes\t{len(gcf_order)}\n")
        fh.write(f"unique_hashes\t{records_written}\n")
        fh.write(f"hash_collisions_across_gcf\t{collisions}\n")
        fh.write(f"length_issues\t{length_issues}\n")
        fh.write(f"enzyme\t{enzyme}\n")

    sys.stderr.write(
        f"Wrote {output_path}: {records_written:,} unique hashes, "
        f"{collisions} cross-GCF collisions\n"
    )


def main():
    ap = argparse.ArgumentParser(
        description="Convert 2bRAD-M pre-digested tag FASTA to fast2bRAD-M .iibdb"
    )
    ap.add_argument("-i", "--input", required=True, help="Input FASTA (.fa/.fa.gz)")
    ap.add_argument("-o", "--output", required=True, help="Output .iibdb path")
    ap.add_argument("-s", "--enzyme", default="BcgI", help="Enzyme name (recorded in stats)")
    ap.add_argument(
        "--expected-length",
        type=int,
        default=None,
        help="Expected tag length; warn if mismatched",
    )
    args = ap.parse_args()

    convert_fasta_to_iibdb(
        Path(args.input),
        Path(args.output),
        args.enzyme,
        args.expected_length,
    )


if __name__ == "__main__":
    main()
