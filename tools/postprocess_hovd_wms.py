#!/usr/bin/env python3
"""
postprocess_hovd_wms.py

在 WMS 8 酶 array job 全部跑完后，收集成功的 per-sample VIP2B profile，
做 merge、HOVD 注释、基础统计，并输出失败样品清单。

用法（HPC 上）:
  python Fast2bRAD-M/tools/postprocess_hovd_wms.py \
    -r results_wms_parallel \
    -d hovd_db \
    -o results_wms_parallel/final \
    -m metadata.tsv.gz
"""
import argparse
import gzip
import os
import shutil
import subprocess
import sys
from pathlib import Path
from typing import Dict, List, Tuple


def find_profiles(results_dir: Path) -> List[Tuple[str, Path]]:
    profiles = []
    for p in sorted(results_dir.rglob("*.VIP2B.xls")):
        sample = p.stem.replace(".VIP2B", "")
        profiles.append((sample, p))
    return profiles


def run_merge(fast2brad_m: str, merge_list: Path, outdir: Path, prefix: str):
    outdir.mkdir(parents=True, exist_ok=True)
    cmd = [
        fast2brad_m,
        "merge",
        "-l", str(merge_list),
        "-o", str(outdir),
        "-p", prefix,
    ]
    sys.stderr.write(f"[merge] {prefix}\n")
    subprocess.run(cmd, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)


def count_detected_taxa(all_xls: Path) -> Tuple[int, int]:
    opener = gzip.open if str(all_xls).endswith(".gz") else open
    n_taxa = 0
    n_nonzero_samples = 0
    with opener(all_xls, "rt", encoding="utf-8") as fh:
        header = fh.readline().rstrip("\n\r").split("\t")
        sample_cols = [h for h in header if h.lstrip("#") not in (
            "Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species", "Strain"
        )]
        for line in fh:
            parts = line.rstrip("\n\r").split("\t")
            if len(parts) < len(header):
                continue
            n_taxa += 1
            for j, col in enumerate(sample_cols):
                idx = header.index(col)
                try:
                    if float(parts[idx]) > 0:
                        n_nonzero_samples += 1
                        break
                except ValueError:
                    pass
    return n_taxa, n_nonzero_samples


def write_summary(
    final_dir: Path,
    profiles: List[Tuple[str, Path]],
    all_samples: List[str],
    failed_samples: List[str],
    all_xls: Path,
):
    n_taxa, n_nonzero = count_detected_taxa(all_xls)
    summary_path = final_dir / "summary.tsv"
    with open(summary_path, "w") as fh:
        fh.write(f"total_samples_in_list\t{len(all_samples)}\n")
        fh.write(f"successful_profiles\t{len(profiles)}\n")
        fh.write(f"failed_samples\t{len(failed_samples)}\n")
        fh.write(f"detected_viral_taxa\t{n_taxa}\n")
        fh.write(f"samples_with_any_detection\t{n_nonzero}\n")
        fh.write("successful\t" + ",".join(sorted(s for s, _ in profiles)) + "\n")
        fh.write("failed\t" + ",".join(sorted(failed_samples)) + "\n")
    sys.stderr.write(f"Wrote {summary_path}\n")


def main():
    ap = argparse.ArgumentParser(description="Post-process HOVD WMS 8-enzyme array results")
    ap.add_argument("-r", "--results-dir", required=True, help="Directory with per-sample results")
    ap.add_argument("-d", "--database", required=True, help="HOVD database directory")
    ap.add_argument("-m", "--metadata", required=True, help="HOVD metadata.tsv.gz")
    ap.add_argument("-o", "--outdir", required=True, help="Final output directory")
    ap.add_argument("-s", "--samples", required=True, help="Original samples_wms.tsv")
    ap.add_argument("--fast2brad-m", default=None, help="Path to fast2bRAD-M binary")
    ap.add_argument("--prefix", default="HOVD_WMS", help="Prefix for merged table")
    ap.add_argument("--annotate-script", default=None, help="Path to annotate_hovd.py")
    args = ap.parse_args()

    results_dir = Path(args.results_dir).resolve()
    outdir = Path(args.outdir).resolve()
    outdir.mkdir(parents=True, exist_ok=True)

    fast2brad_m = args.fast2brad_m
    if not fast2brad_m:
        repo_bin = Path(__file__).resolve().parent.parent / "target" / "release" / "fast2bRAD-M"
        fast2brad_m = shutil.which("fast2bRAD-M") or str(repo_bin)
    if not Path(fast2brad_m).exists():
        sys.stderr.write(f"ERROR: fast2bRAD-M not found: {fast2brad_m}\n")
        sys.exit(1)

    # All samples in original list
    all_samples = []
    with open(args.samples, "rt", encoding="utf-8") as fh:
        for line in fh:
            line = line.rstrip("\n\r").strip()
            if not line or line.startswith("#"):
                continue
            all_samples.append(line.split("\t")[0].strip())

    profiles = find_profiles(results_dir)
    successful = {s for s, _ in profiles}
    failed_samples = [s for s in all_samples if s not in successful]

    if not profiles:
        sys.stderr.write("ERROR: no successful VIP2B profiles found\n")
        sys.exit(1)

    merge_list = outdir / "merge.list"
    with open(merge_list, "w") as fh:
        for sample, path in profiles:
            fh.write(f"{sample}\t{path}\n")

    run_merge(fast2brad_m, merge_list, outdir, args.prefix)
    all_xls = outdir / f"{args.prefix}.all.xls"
    filtered_xls = outdir / f"{args.prefix}.filtered.xls"

    annotate_script = args.annotate_script
    if not annotate_script:
        annotate_script = Path(__file__).resolve().parent / "annotate_hovd.py"
    annotate_cmd = [
        sys.executable,
        str(annotate_script),
        "-i", str(all_xls),
        "-d", str(args.metadata),
        "-o", str(outdir / "annotation"),
    ]
    sys.stderr.write("[annotate] running annotate_hovd.py\n")
    subprocess.run(annotate_cmd, check=True)

    write_summary(outdir, profiles, all_samples, failed_samples, all_xls)

    # Also write failed list for easy downstream handling
    with open(outdir / "failed_samples.txt", "w") as fh:
        fh.write("\n".join(failed_samples) + "\n")

    sys.stderr.write("Post-processing complete.\n")


if __name__ == "__main__":
    main()
