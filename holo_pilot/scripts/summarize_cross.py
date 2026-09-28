#!/usr/bin/env python3
"""Summarize human-vs-microbial tag collision results from f2brad-host cross.

Inputs (per enzyme):
  cross_results/{enzyme}/
    collisions.{enzyme}.2.tsv
    species_summary.{enzyme}.2.tsv
    human_mask.{enzyme}.2.txt

Output:
  cross_collision_summary.md
"""

import os
import sys
from collections import defaultdict
from pathlib import Path

ENZYMES = ["BcgI", "BsaXI", "AlfI"]
UNIQUE_TAGS = {"BcgI": 108697, "BsaXI": 790007, "AlfI": 359400}
RESULTS_DIR = Path(__file__).parent.parent / "cross_results"
OUT_MD = Path(__file__).parent.parent / "results" / "cross_collision_summary.md"


def read_mask_count(path: Path) -> int:
    if not path.exists():
        return 0
    with open(path) as f:
        return sum(1 for line in f if line.strip())


def read_cross_summary(path: Path) -> dict:
    """Parse cross_summary.{enzyme}.2.txt into a dict."""
    summary = {}
    if not path.exists():
        return summary
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line or "\t" not in line:
                continue
            key, value = line.split("\t", 1)
            try:
                summary[key] = int(value)
            except ValueError:
                try:
                    summary[key] = float(value)
                except ValueError:
                    summary[key] = value
    return summary


def read_species_summary(path: Path):
    """Return list of (species, c0, c1, c2, total) sorted by total descending."""
    rows = []
    if not path.exists():
        return rows
    with open(path) as f:
        next(f)  # header
        for line in f:
            line = line.strip()
            if not line:
                continue
            parts = line.split("\t")
            if len(parts) < 5:
                continue
            species, c0, c1, c2, total = parts[:5]
            rows.append((species, int(c0), int(c1), int(c2), int(total)))
    rows.sort(key=lambda x: x[4], reverse=True)
    return rows


def read_collision_events(path: Path) -> dict:
    """Return dict distance -> count of collision events."""
    counts = defaultdict(int)
    if not path.exists():
        return counts
    with open(path) as f:
        next(f)  # header
        for line in f:
            line = line.strip()
            if not line:
                continue
            parts = line.split("\t")
            if len(parts) < 3:
                continue
            dist = int(parts[2])
            counts[dist] += 1
    return counts


def main():
    sections = []
    sections.append("# Human↔microbial tag cross-assignment collision summary\n")
    sections.append(
        "Reference: T2T-CHM13v2.0 human tags vs. GTDB representative genomes "
        "(one per species, 143,586 genomes).\n\n"
    )
    sections.append("Tool: `f2brad-host cross` (Fast2bRAD-M `holo` branch)\n\n")
    sections.append("Max Hamming distance for collision: 2\n\n")

    # Table 1: mask fraction per enzyme
    sections.append("## Human tag masking fraction\n")
    sections.append(
        "A human unique tag is masked if any microbial tag lies within ≤2 "
        "mismatches of it. Masked tags should be excluded from the host "
        "genotyping panel to prevent microbial reads from being mis-assigned "
        "as host alleles.\n\n"
    )
    sections.append("| enzyme | unique human tags | masked tags | mask fraction |\n")
    sections.append("|--------|------------------:|------------:|--------------:|\n")

    microbe_table = ["## Microbial tag masking fraction\n"]
    microbe_table.append(
        "A microbial tag is masked if it lies within ≤2 mismatches of any human "
        "tag. Masked tags should be excluded from the microbial reference "
        "database so that host reads are not mis-assigned as microbial.\n\n"
    )
    microbe_table.append(
        "| enzyme | total microbial tags | unique microbial tags | masked microbial tags | microbial mask fraction |\n"
    )
    microbe_table.append(
        "|--------|---------------------:|----------------------:|----------------------:|------------------------:|\n"
    )

    mismatch_table = ["## Mismatch-distance breakdown\n"]
    mismatch_table.append(
        "Collision events are counted by Hamming distance. "
        "0-mismatch = identical tag sequence; 1/2-mismatch = near-identical.\n\n"
    )
    mismatch_table.append("| enzyme | 0-mismatch | 1-mismatch | 2-mismatch | total events |\n")
    mismatch_table.append("|--------|-----------:|-----------:|-----------:|-------------:|\n")

    species_blocks = ["## Top colliding microbial species\n"]
    species_blocks.append(
        "Species are ranked by total collision events (sum of 0/1/2 mismatches).\n\n"
    )

    for enzyme in ENZYMES:
        base = RESULTS_DIR / enzyme
        mask_path = base / f"human_mask.{enzyme}.2.txt"
        microbe_mask_path = base / f"microbe_mask.{enzyme}.2.txt"
        summary_path = base / f"species_summary.{enzyme}.2.tsv"
        collisions_path = base / f"collisions.{enzyme}.2.tsv"
        cross_summary_path = base / f"cross_summary.{enzyme}.2.txt"

        unique = UNIQUE_TAGS[enzyme]
        masked = read_mask_count(mask_path)
        frac = masked / unique if unique else 0.0
        sections.append(
            f"| {enzyme} | {unique:,} | {masked:,} | {frac:.4%} |\n"
        )

        cross_summary = read_cross_summary(cross_summary_path)
        total_microbe = cross_summary.get("total_microbial_tags", 0)
        unique_microbe = cross_summary.get("unique_microbial_tags", 0)
        masked_microbe = read_mask_count(microbe_mask_path)
        microbe_frac = (
            masked_microbe / unique_microbe if unique_microbe else 0.0
        )
        microbe_table.append(
            f"| {enzyme} | {total_microbe:,} | {unique_microbe:,} | {masked_microbe:,} | {microbe_frac:.4%} |\n"
        )

        event_counts = read_collision_events(collisions_path)
        c0 = event_counts.get(0, 0)
        c1 = event_counts.get(1, 0)
        c2 = event_counts.get(2, 0)
        total = c0 + c1 + c2
        mismatch_table.append(
            f"| {enzyme} | {c0:,} | {c1:,} | {c2:,} | {total:,} |\n"
        )

        species_rows = read_species_summary(summary_path)
        species_blocks.append(f"### {enzyme}\n")
        species_blocks.append("| rank | species | 0-mismatch | 1-mismatch | 2-mismatch | total |\n")
        species_blocks.append("|-----:|---------|-----------:|-----------:|-----------:|------:|\n")
        for i, (sp, s0, s1, s2, st) in enumerate(species_rows[:10], 1):
            species_blocks.append(
                f"| {i} | {sp} | {s0:,} | {s1:,} | {s2:,} | {st:,} |\n"
            )
        if not species_rows:
            species_blocks.append("_No collisions found._\n")
        species_blocks.append("\n")

    sections.append("\n")
    sections.extend(microbe_table)
    sections.append("\n")
    sections.extend(mismatch_table)
    sections.append("\n")
    sections.extend(species_blocks)

    # Interpretation — updated once real microbial-tag counts are available.
    sections.append("## Interpretation\n\n")
    sections.append(
        "- **Microbial mask fractions are negligible.** Even though GTDB representative "
        "genomes contribute hundreds of millions of unique tags per enzyme, only a "
        "tiny fraction collide with human tags: **0.0002% for BcgI**, **0.0003% for "
        "AlfI**, and **0.0517% for BsaXI**. Removing these tags from the microbial "
        "reference database costs almost nothing in panel size.\n"
    )
    sections.append(
        "- **Total microbial tag counts differ sharply from human.** In the human genome "
        "BcgI is CpG-penalised and yields the fewest sites, but microbes have no CpG "
        "depletion, so BcgI produces the **largest microbial tag pool (~386 M total, "
        "~332 M unique)**. BsaXI and AlfI give smaller microbial pools (~190 M and "
        "~167 M total, respectively) because their motifs are less common in microbes.\n"
    )
    sections.append(
        "- **BsaXI is the collision outlier in both directions.** It has both the highest "
        "human mask fraction (7.84%) and the highest microbial mask fraction "
        "(~0.052%). The microbial-arm cost is still small, but the host-arm cost is "
        "material and must be reported as a limitation of the shorter 27 bp tag.\n"
    )
    sections.append(
        "- **0-mismatch collisions are the most dangerous:** they are indistinguishable "
        "from true host tags at the sequence level. 1/2-mismatch collisions are still "
        "detectable if the matching engine reports distance, but they inflate genotype "
        "uncertainty.\n"
    )
    sections.append(
        "- **Both masks are required for a holo pipeline.** The human mask protects host "
        "genotyping from microbial reads; the microbial mask protects abundance "
        "estimation from host reads. Because the masks are built from the same cross "
        "collision list, they are mutually consistent.\n"
    )

    OUT_MD.parent.mkdir(parents=True, exist_ok=True)
    with open(OUT_MD, "w") as f:
        f.writelines(sections)
    print(f"Wrote {OUT_MD}")


if __name__ == "__main__":
    main()
