#!/usr/bin/env python3
"""Build microbiome + host feature tables from f2brad-holo diurnal outputs."""
import argparse
from pathlib import Path

import numpy as np
import pandas as pd

CANDIDATE_INTERVALS = {
    "FUT2": ("NC_060943.1", 51690377, 51700359),
    "LCT": ("NC_060926.1", 136228325, 136281636),
    "ABO": ("NC_060933.1", 145463984, 145489076),
}


def load_species_counts(holo_dir: Path, sample: str) -> pd.Series:
    path = holo_dir / sample / "species_counts.tsv"
    if not path.exists():
        return pd.Series(dtype=float)
    df = pd.read_csv(path, sep="\t")
    if df.empty:
        return pd.Series(dtype=float)
    return df.set_index("species")["count"]


def load_host_fraction(holo_dir: Path, sample: str) -> float:
    path = holo_dir / sample / "holo_classify.tsv"
    if not path.exists():
        return float("nan")
    df = pd.read_csv(path, sep="\t")
    return df.set_index("metric").loc["host_fraction", "value"]


def parse_candidate_dosages(holo_dir: Path, sample: str):
    path = holo_dir / sample / "genotypes.vcf"
    if not path.exists():
        return {}
    dosages = {}
    with open(path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            parts = line.split("\t")
            chrom, pos = parts[0], int(parts[1])
            fmt, sample_field = parts[8], parts[9]
            fmt_map = dict(zip(fmt.split(":"), sample_field.split(":")))
            gt = fmt_map.get("GT", "./.")
            dosage = -1.0
            if gt == "0/0":
                dosage = 0.0
            elif gt == "0/1":
                dosage = 1.0
            elif gt == "1/1":
                dosage = 2.0
            for gene, (gchrom, start, end) in CANDIDATE_INTERVALS.items():
                if chrom == gchrom and start <= pos <= end:
                    dosages[f"{gene}_{chrom}_{pos}"] = dosage
    return dosages


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--metadata", required=True)
    parser.add_argument("--holo-dir", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    outdir = Path(args.output)
    outdir.mkdir(parents=True, exist_ok=True)

    meta = pd.read_csv(args.metadata, sep="\t")
    holo_root = Path(args.holo_dir)

    species_list = []
    host_rows = []
    for _, row in meta.iterrows():
        sample = row["sample_id"]
        counts = load_species_counts(holo_root, sample)
        species_list.append(counts.rename(sample))

        host_frac = load_host_fraction(holo_root, sample)
        snp_dosages = parse_candidate_dosages(holo_root, sample)
        host_row = {"sample_id": sample, "host_fraction": host_frac}
        host_row.update(snp_dosages)
        host_rows.append(host_row)

    species_df = pd.concat(species_list, axis=1).fillna(0).T
    # CLR-like transform: add pseudo-count, log
    species_clr = np.log1p(species_df.div(species_df.sum(axis=1), axis=0) * 1e6)

    host_df = pd.DataFrame(host_rows).set_index("sample_id").fillna(-1)

    species_clr = species_clr.loc[meta["sample_id"]]
    host_df = host_df.loc[meta["sample_id"]]

    species_clr.to_csv(outdir / "microbiome_clr.tsv", sep="\t")
    host_df.to_csv(outdir / "host_features.tsv", sep="\t")
    meta.set_index("sample_id").to_csv(outdir / "sample_metadata.tsv", sep="\t")

    print(f"Wrote {species_clr.shape[1]} microbiome features and "
          f"{host_df.shape[1]} host features for {len(meta)} samples")


if __name__ == "__main__":
    main()
