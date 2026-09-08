#!/usr/bin/env python3
"""
quantify_vip2b.py

使用 Fast2bRAD-M 的 extract/quantify 子命令，对 VIP2B 风格的 8 酶组合
病毒数据库进行样品定量。

流程:
  1. 对每对 (R1,R2) 或 SE reads，分别用 8 种酶做 fast2bRAD-M extract，得到
     8 个 .iibsp 文件。
  2. 将 8 个 .iibsp 合并为一个 combined .iibsp（病毒 tag 集合）。
  3. 用 fast2bRAD-M quantify + 转换后的 VIP2B 数据库得到每个样品的
     病毒丰度表。
  4. （可选）用 fast2bRAD-M merge 输出跨样品的丰度矩阵。

输入样品列表格式（tab 分隔）:
  sample_name  R1.fq.gz  [R2.fq.gz]

数据库目录要求（由 convert_vip2b_db.py 生成）:
  {db_dir}/{site}.{level}.iibdb
  {db_dir}/abfh_classify_with_speciename.txt.gz

用法示例:
  python tools/quantify_vip2b.py \
    -i samples.tsv \
    -d vip2b_db/ \
    -o vip2b_results/ \
    -j 16
"""
import argparse
import gzip
import os
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path
from typing import List, Optional, Tuple

DEFAULT_ENZYMES = ["AlfI", "BcgI", "BslFI", "CjeI", "CjePI", "FalI", "HaeIV", "Hin4I"]
SAMPLE_TAG_MAGIC = b"IIBS"
SAMPLE_TAG_VERSION = 1


def find_fast2brad_m(args_bin: Optional[str]) -> str:
    if args_bin:
        return os.path.abspath(args_bin)
    # 优先使用仓库 release 二进制，否则在 PATH 中查找
    repo_bin = Path(__file__).resolve().parent.parent / "target" / "release" / "fast2bRAD-M"
    if repo_bin.exists():
        return str(repo_bin)
    path_bin = shutil.which("fast2bRAD-M")
    if path_bin:
        return path_bin
    sys.stderr.write("ERROR: fast2bRAD-M binary not found. Use --fast2brad-m.\n")
    sys.exit(1)


def parse_sample_list(path: Path) -> List[Tuple[str, Path, Optional[Path]]]:
    samples = []
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt", encoding="utf-8") as fh:
        for line in fh:
            line = line.rstrip("\n\r").strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t")
            name = parts[0].strip()
            r1 = Path(parts[1].strip())
            r2 = Path(parts[2].strip()) if len(parts) > 2 and parts[2].strip() else None
            if not r1.exists():
                sys.stderr.write(f"Warning: R1 not found for {name}: {r1}\n")
                continue
            if r2 and not r2.exists():
                sys.stderr.write(f"Warning: R2 not found for {name}: {r2}\n")
                continue
            samples.append((name, r1, r2))
    return samples


def run_extract(
    fast2brad_m: str,
    enzyme: str,
    r1: Path,
    r2: Optional[Path],
    outdir: Path,
    sample: str,
    input_type: int,
    qc: str,
    qc_scope: str,
    threads: int,
):
    outdir.mkdir(parents=True, exist_ok=True)
    cmd = [
        fast2brad_m,
        "extract",
        "-t", str(input_type),
        "-s", enzyme,
        "--od", str(outdir),
        "--op", sample,
        "-j", str(threads),
        "--qc", qc,
        "--qc-scope", qc_scope,
        "-i", str(r1),
    ]
    if r2:
        cmd.append(str(r2))
    sys.stderr.write(f"[extract] {sample} {enzyme}\n")
    subprocess.run(cmd, check=True)


def merge_iibsp(inputs: List[Path], output: Path):
    """将多个 .iibsp 文件合并成一个（跳过每个文件的 8 字节 header）。"""
    output.parent.mkdir(parents=True, exist_ok=True)
    with open(output, "wb") as outfh:
        outfh.write(SAMPLE_TAG_MAGIC)
        outfh.write(SAMPLE_TAG_VERSION.to_bytes(4, "little"))
        for p in inputs:
            with open(p, "rb") as infh:
                header = infh.read(8)
                if len(header) != 8 or header[:4] != SAMPLE_TAG_MAGIC:
                    sys.stderr.write(f"Warning: {p} does not look like an .iibsp file; copying raw\n")
                    if len(header) == 8:
                        outfh.write(header)
                shutil.copyfileobj(infh, outfh)


def quantify_sample(
    fast2brad_m: str,
    sample: str,
    combined_iibsp: Path,
    db_dir: Path,
    level: str,
    site: str,
    outdir: Path,
    threads: int,
    gscore: float,
):
    list_file = outdir / "samples.iibsp.list"
    list_file.parent.mkdir(parents=True, exist_ok=True)
    with open(list_file, "w") as fh:
        fh.write(f"{sample}\t{combined_iibsp}\n")

    cmd = [
        fast2brad_m,
        "quantify",
        "-l", str(list_file),
        "-d", str(db_dir),
        "-t", level,
        "-s", site,
        "-o", str(outdir),
        "-j", str(threads),
        "-g", str(gscore),
        "-v", "yes",
    ]
    sys.stderr.write(f"[quantify] {sample}\n")
    subprocess.run(cmd, check=True)
    # quantify 输出在 outdir/{sample}/{sample}.{site}.xls
    return outdir / sample / f"{sample}.{site}.xls"


def run_merge(
    fast2brad_m: str,
    merge_list: Path,
    outdir: Path,
    prefix: str,
):
    outdir.mkdir(parents=True, exist_ok=True)
    cmd = [
        fast2brad_m,
        "merge",
        "-l", str(merge_list),
        "-o", str(outdir),
        "-p", prefix,
    ]
    sys.stderr.write("[merge] merging per-sample profiles\n")
    subprocess.run(cmd, check=True)


def main():
    ap = argparse.ArgumentParser(description="VIP2B-style 8-enzyme viral profiling with Fast2bRAD-M")
    ap.add_argument("-i", "--input", required=True, help="Sample list TSV: sample<TAB>R1<TAB>[R2]")
    ap.add_argument("-d", "--database", required=True, help="VIP2B database directory")
    ap.add_argument("-o", "--outdir", required=True, help="Output directory")
    ap.add_argument("-l", "--level", default="species", help="Taxonomy level (species/genus/...)")
    ap.add_argument("-s", "--site", default="BcgI", help="Placeholder enzyme/site name used in convert_vip2b_db.py")
    ap.add_argument("-e", "--enzymes", default=",".join(DEFAULT_ENZYMES), help="Comma-separated enzyme list")
    ap.add_argument("-j", "--threads", type=int, default=4, help="Threads passed to fast2bRAD-M")
    ap.add_argument("--fast2brad-m", default=None, help="Path to fast2bRAD-M binary")
    ap.add_argument("--input-type", type=int, default=2, choices=[2, 3], help="fast2bRAD-M extract input type (2=shotgun, 3=single 2bRAD tags)")
    ap.add_argument("--qc", default="no", help="fast2bRAD-M extract QC (yes/no)")
    ap.add_argument("--qc-scope", default="auto", help="fast2bRAD-M extract QC scope (read/tag/auto)")
    ap.add_argument("-g", "--gscore", type=float, default=0.0, help="G-score threshold for quantify")
    ap.add_argument("--no-merge", action="store_true", help="Skip the final merge step")
    ap.add_argument("--merge-prefix", default="VIP2B", help="Prefix for merged abundance table")
    args = ap.parse_args()

    fast2brad_m = find_fast2brad_m(args.fast2brad_m)
    sys.stderr.write(f"Using fast2bRAD-M: {fast2brad_m}\n")

    samples = parse_sample_list(Path(args.input))
    if not samples:
        sys.stderr.write("ERROR: no valid samples found\n")
        sys.exit(1)

    db_dir = Path(args.database).resolve()
    outdir = Path(args.outdir).resolve()
    outdir.mkdir(parents=True, exist_ok=True)
    enzymes = [e.strip() for e in args.enzymes.split(",") if e.strip()]

    db_file = db_dir / f"{args.site}.{args.level}.iibdb"
    classify_file = db_dir / "abfh_classify_with_speciename.txt.gz"
    if not db_file.exists():
        sys.stderr.write(f"ERROR: database file not found: {db_file}\n")
        sys.exit(1)
    if not classify_file.exists():
        sys.stderr.write(f"ERROR: classify file not found: {classify_file}\n")
        sys.exit(1)

    extract_root = outdir / "01_extract"
    quantify_root = outdir / "02_quantify"
    profiles_dir = outdir / "03_profiles"
    profiles_dir.mkdir(parents=True, exist_ok=True)

    merge_entries = []
    for sample, r1, r2 in samples:
        sample_extract_dir = extract_root / sample
        per_enzyme_iibsp = []
        for enzyme in enzymes:
            run_extract(
                fast2brad_m,
                enzyme,
                r1,
                r2,
                sample_extract_dir,
                sample,
                args.input_type,
                args.qc,
                args.qc_scope,
                args.threads,
            )
            per_enzyme_iibsp.append(sample_extract_dir / f"{sample}.{enzyme}.iibsp")
            # 检查文件是否存在
            if not per_enzyme_iibsp[-1].exists():
                sys.stderr.write(f"ERROR: extract failed for {sample} {enzyme}\n")
                sys.exit(1)

        combined = sample_extract_dir / f"{sample}.VIP2B.iibsp"
        merge_iibsp(per_enzyme_iibsp, combined)

        xls_path = quantify_sample(
            fast2brad_m,
            sample,
            combined,
            db_dir,
            args.level,
            args.site,
            quantify_root,
            args.threads,
            args.gscore,
        )
        if not xls_path.exists():
            sys.stderr.write(
                f"Warning: no abundance table produced for {sample} (no viral tags detected)\n"
            )
            continue
        # 创建易访问的拷贝
        profile_copy = profiles_dir / f"{sample}.VIP2B.xls"
        shutil.copy2(xls_path, profile_copy)
        merge_entries.append((sample, profile_copy))

    if not args.no_merge and len(merge_entries) > 1:
        merge_list = outdir / "merge.list"
        with open(merge_list, "w") as fh:
            for sample, path in merge_entries:
                fh.write(f"{sample}\t{path}\n")
        run_merge(fast2brad_m, merge_list, outdir, args.merge_prefix)
        sys.stderr.write(f"Merged tables written to {outdir}/{args.merge_prefix}.{{all,filtered}}.xls\n")

    sys.stderr.write("All done.\n")


if __name__ == "__main__":
    main()
