#!/usr/bin/env python3
"""
quantify_vip2b_parallel.py

与 quantify_vip2b.py 相同功能，但把 8 个酶的 extract 在同一节点上并行运行，
避免顺序读取 8 次大 WMS 文件造成的总时间过长。

用法:
  python tools/quantify_vip2b_parallel.py \
    -i samples.tsv \
    -d hovd_db/ \
    -o hovd_results/ \
    -j 16
"""
import argparse
import gzip
import os
import shutil
import subprocess
import sys
import tempfile
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
from typing import List, Optional, Tuple

DEFAULT_ENZYMES = ["AlfI", "BcgI", "BslFI", "CjeI", "CjePI", "FalI", "HaeIV", "Hin4I"]
SAMPLE_TAG_MAGIC = b"IIBS"
SAMPLE_TAG_VERSION = 1


def find_fast2brad_m(args_bin: Optional[str]) -> str:
    if args_bin:
        return os.path.abspath(args_bin)
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
) -> Path:
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
    subprocess.run(cmd, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    return outdir / f"{sample}.{enzyme}.iibsp"


def merge_iibsp(inputs: List[Path], output: Path):
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
    subprocess.run(cmd, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    return outdir / sample / f"{sample}.{site}.xls"


def process_one_sample(
    fast2brad_m: str,
    db_dir: Path,
    sample: str,
    r1: Path,
    r2: Optional[Path],
    outdir: Path,
    enzymes: List[str],
    input_type: int,
    qc: str,
    qc_scope: str,
    extract_threads: int,
    quant_threads: int,
    gscore: float,
    site: str,
    level: str,
) -> Optional[Path]:
    sample_extract_dir = outdir / "01_extract" / sample / sample
    sample_extract_dir.mkdir(parents=True, exist_ok=True)

    # Run 8 enzyme extracts in parallel (use threads because work is subprocess-based)
    per_enzyme_iibsp = []
    with ThreadPoolExecutor(max_workers=len(enzymes)) as executor:
        futures = {
            executor.submit(
                run_extract,
                fast2brad_m,
                enzyme,
                r1,
                r2,
                sample_extract_dir,
                sample,
                input_type,
                qc,
                qc_scope,
                extract_threads,
            ): enzyme
            for enzyme in enzymes
        }
        for future in as_completed(futures):
            enzyme = futures[future]
            try:
                path = future.result()
                per_enzyme_iibsp.append(path)
            except Exception as e:
                sys.stderr.write(f"ERROR: extract failed for {sample} {enzyme}: {e}\n")
                return None

    if len(per_enzyme_iibsp) != len(enzymes):
        sys.stderr.write(f"ERROR: not all enzymes succeeded for {sample}\n")
        return None

    combined = sample_extract_dir.parent / f"{sample}.VIP2B.iibsp"
    merge_iibsp(per_enzyme_iibsp, combined)

    quantify_root = outdir / "02_quantify"
    xls_path = quantify_sample(
        fast2brad_m,
        sample,
        combined,
        db_dir,
        level,
        site,
        quantify_root,
        quant_threads,
        gscore,
    )
    if not xls_path.exists():
        sys.stderr.write(f"Warning: no abundance table produced for {sample}\n")
        return None

    profiles_dir = outdir / "03_profiles"
    profiles_dir.mkdir(parents=True, exist_ok=True)
    profile_copy = profiles_dir / f"{sample}.VIP2B.xls"
    shutil.copy2(xls_path, profile_copy)
    return profile_copy


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
    ap = argparse.ArgumentParser(description="VIP2B-style 8-enzyme viral profiling with parallel enzyme extraction")
    ap.add_argument("-i", "--input", required=True, help="Sample list TSV: sample<TAB>R1<TAB>[R2]")
    ap.add_argument("-d", "--database", required=True, help="VIP2B database directory")
    ap.add_argument("-o", "--outdir", required=True, help="Output directory")
    ap.add_argument("-l", "--level", default="species", help="Taxonomy level")
    ap.add_argument("-s", "--site", default="BcgI", help="Placeholder enzyme/site name")
    ap.add_argument("-e", "--enzymes", default=",".join(DEFAULT_ENZYMES), help="Comma-separated enzyme list")
    ap.add_argument("-j", "--threads", type=int, default=16, help="Total threads per sample (split across enzymes)")
    ap.add_argument("--fast2brad-m", default=None, help="Path to fast2bRAD-M binary")
    ap.add_argument("--input-type", type=int, default=2, choices=[2, 3], help="fast2bRAD-M extract input type")
    ap.add_argument("--qc", default="no", help="fast2bRAD-M extract QC")
    ap.add_argument("--qc-scope", default="auto", help="fast2bRAD-M extract QC scope")
    ap.add_argument("-g", "--gscore", type=float, default=0.0, help="G-score threshold")
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

    # Each enzyme gets 2 threads (or 1 if total < enzymes)
    extract_threads = max(1, args.threads // len(enzymes))
    quant_threads = max(4, args.threads)

    merge_entries = []
    failed_samples: list[str] = []
    for sample, r1, r2 in samples:
        profile = process_one_sample(
            fast2brad_m,
            db_dir,
            sample,
            r1,
            r2,
            outdir,
            enzymes,
            args.input_type,
            args.qc,
            args.qc_scope,
            extract_threads,
            quant_threads,
            args.gscore,
            args.site,
            args.level,
        )
        if profile:
            merge_entries.append((sample, profile))
        else:
            failed_samples.append(sample)

    if not args.no_merge and len(merge_entries) > 1:
        merge_list = outdir / "merge.list"
        with open(merge_list, "w") as fh:
            for sample, path in merge_entries:
                fh.write(f"{sample}\t{path}\n")
        run_merge(fast2brad_m, merge_list, outdir, args.merge_prefix)
        sys.stderr.write(f"Merged tables written to {outdir}/{args.merge_prefix}.{{all,filtered}}.xls\n")

    if failed_samples:
        sys.stderr.write(f"ERROR: {len(failed_samples)} sample(s) failed: {', '.join(failed_samples)}\n")
        sys.exit(1)

    sys.stderr.write("All done.\n")


if __name__ == "__main__":
    main()
