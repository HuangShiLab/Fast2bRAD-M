#!/usr/bin/env python3
"""Build microbiome + host feature tables from f2brad-holo ECC outputs."""
import argparse
import csv
import re
from pathlib import Path

import pandas as pd

# Candidate gene intervals on T2T-CHM13v2.0 RefSeq (1-based inclusive).
# Includes classic caries-susceptibility genes (enamel, immune, taste, vitamin D)
# plus milk-digestion loci as positive controls.
CANDIDATE_INTERVALS = {
    # Positive controls: diet/lactase related
    "FUT2": ("NC_060943.1", 51690377, 51700359),
    "LCT": ("NC_060926.1", 136228325, 136281636),
    "ABO": ("NC_060933.1", 145463984, 145489076),
    # Enamel / dentin matrix genes
    "ENAM": ("NC_060928.1", 73971213, 73989293),
    "AMBN": ("NC_060928.1", 73933220, 73948252),
    "MMP20": ("NC_060935.1", 102579028, 102627492),
    "AMELX": ("NC_060947.1", 10875912, 10892087),
    "TUFT1": ("NC_060925.1", 150663979, 150707247),
    "KLK4": ("NC_060943.1", 53994435, 53999407),
    "DSPP": ("NC_060928.1", 90935096, 90943283),
    # Immune / defence / taste / vitamin D
    "TAS2R38": ("NC_060931.1", 143288364, 143289506),
    "DEFB1": ("NC_060932.1", 6625781, 6633102),
    "MBL2": ("NC_060934.1", 53612365, 53619761),
    "VDR": ("NC_060936.1", 47802918, 47866398),
    # Matrix metalloproteinases
    "MMP2": ("NC_060940.1", 61276956, 61304813),
    "MMP3": ("NC_060935.1", 102839570, 102847378),
    "MMP9": ("NC_060944.1", 47744921, 47752589),
    "MMP13": ("NC_060935.1", 102946765, 102959501),
}


def load_metadata(path: str) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t")
    required = {"sample_id", "r1", "r2", "phenotype"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"metadata missing columns: {missing}")
    return df


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


def parse_vcf_for_candidates(holo_dir: Path, sample: str, gene_intervals: dict):
    """Return a dict of candidate-gene SNP dosages (0/1/2 or -1 for missing)."""
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
            for gene, (gchrom, start, end) in gene_intervals.items():
                if chrom == gchrom and start <= pos <= end:
                    dosages[f"{gene}_{chrom}_{pos}"] = dosage
    return dosages


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--metadata", required=True)
    parser.add_argument("--holo-dir", required=True)
    parser.add_argument("--enzyme", default="BcgI")
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    outdir = Path(args.output)
    outdir.mkdir(parents=True, exist_ok=True)

    meta = load_metadata(args.metadata)
    holo_root = Path(args.holo_dir)

    # Microbiome features
    species_list = []
    host_rows = []
    for _, row in meta.iterrows():
        sample = row["sample_id"]
        counts = load_species_counts(holo_root, sample)
        species_list.append(counts.rename(sample))

        host_frac = load_host_fraction(holo_root, sample)
        snp_dosages = parse_vcf_for_candidates(holo_root, sample, CANDIDATE_INTERVALS)
        host_row = {"sample_id": sample, "host_fraction": host_frac}
        host_row.update(snp_dosages)
        host_rows.append(host_row)

    species_df = pd.concat(species_list, axis=1).fillna(0).T
    # Proportion normalization
    species_df = species_df.div(species_df.sum(axis=1), axis=0).fillna(0)

    host_df = pd.DataFrame(host_rows).set_index("sample_id")
    # fill missing SNPs with dosage -1 for downstream imputation
    host_df = host_df.fillna(-1)

    # Keep only polymorphic SNPs: at least 2 non-missing samples and
    # more than one distinct dosage among non-missing calls.
    snp_cols = [c for c in host_df.columns if c != "host_fraction"]
    keep_cols = ["host_fraction"]
    for col in snp_cols:
        observed = host_df[col].replace(-1, pd.NA).dropna().astype(float)
        if observed.nunique() > 1 and observed.shape[0] >= 2:
            keep_cols.append(col)
    host_df = host_df[keep_cols]

    # Align to metadata order
    species_df = species_df.loc[meta["sample_id"]]
    host_df = host_df.loc[meta["sample_id"]]

    species_df.to_csv(outdir / "X_microbiome.tsv", sep="\t")
    host_df.to_csv(outdir / "X_host.tsv", sep="\t")
    meta[["sample_id", "phenotype"]].set_index("sample_id").to_csv(
        outdir / "y.tsv", sep="\t"
    )

    print(f"Wrote {species_df.shape[1]} microbiome features and "
          f"{host_df.shape[1]} host features for {len(meta)} samples")


if __name__ == "__main__":
    main()
