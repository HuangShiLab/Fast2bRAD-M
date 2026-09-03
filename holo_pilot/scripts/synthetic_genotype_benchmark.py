#!/usr/bin/env python3
"""Synthetic benchmark for f2brad-host genotype accuracy.

Reads the host tag DB, simulates paired-end fragments from a subset of loci
with known hom-ref/het/hom-alt genotypes and random sequencing errors, then
runs f2brad-host genotype and compares the emitted VCF to the truth.
"""
import argparse
import csv
import gzip
import math
import os
import random
import subprocess
import sys
from collections import defaultdict
from pathlib import Path

BASES = [b'A', b'C', b'G', b'T']
COMP = {b'A': b'T', b'C': b'G', b'G': b'C', b'T': b'A',
        b'a': b't', b'c': b'g', b'g': b'c', b't': b'a'}


def reverse_complement(seq: bytes) -> bytes:
    return bytes(COMP.get(b, b) for b in reversed(seq))


def mutate_base(base: bytes) -> bytes:
    return random.choice([b for b in BASES if b != base])


def add_errors(seq: bytes, error_rate: float, error_qual: int) -> tuple[bytes, bytes]:
    """Introduce random errors and return (seq, qual)."""
    seq = bytearray(seq)
    qual = bytearray()
    for i, base in enumerate(seq):
        if random.random() < error_rate:
            seq[i] = mutate_base(bytes([base]))[0]
            qual.append(error_qual + 33)
        else:
            qual.append(40 + 33)  # Phred 40
    return bytes(seq), bytes(qual)


def load_host_db(path: str, n_loci: int):
    loci = []
    opener = gzip.open if path.endswith('.gz') else open
    with opener(path, 'rt') as fh:
        reader = csv.DictReader(fh, delimiter='\t')
        for row in reader:
            loci.append({
                'contig': row['contig'],
                'pos': int(row['pos']),
                'seq': row['seq'].encode(),
                'canonical': row['canonical'].encode(),
            })
            if len(loci) >= n_loci:
                break
    return loci


def simulate(loci, coverage, error_rate, error_qual, seed):
    random.seed(seed)
    truth = {}
    r1_records = []
    r2_records = []
    frag_id = 0
    for locus_idx, locus in enumerate(loci):
        ref_tag = locus['seq']
        gt = random.choice(['0/0', '0/1', '1/1'])
        # choose a single tag position for the alt allele
        alt_pos = len(ref_tag) // 2
        alt_base = mutate_base(bytes([ref_tag[alt_pos]]))
        truth[(locus_idx, alt_pos)] = {
            'contig': locus['contig'],
            'genomic_pos': locus['pos'] + alt_pos,
            'ref': bytes([ref_tag[alt_pos]]),
            'alt': alt_base,
            'gt': gt,
        }
        for _ in range(coverage):
            # draw allele for this fragment
            if gt == '0/0':
                is_alt = False
            elif gt == '1/1':
                is_alt = True
            else:
                is_alt = random.random() < 0.5
            frag_tag = bytearray(ref_tag)
            if is_alt:
                frag_tag[alt_pos] = alt_base[0]
            frag_tag = bytes(frag_tag)

            r1_seq, r1_qual = add_errors(frag_tag, error_rate, error_qual)
            r2_seq, r2_qual = add_errors(reverse_complement(frag_tag), error_rate, error_qual)

            name = f"frag{frag_id}"
            r1_records.append((name, r1_seq, r1_qual))
            r2_records.append((name, r2_seq, r2_qual))
            frag_id += 1
    return truth, r1_records, r2_records


def write_fastq(path: str, records):
    opener = gzip.open if path.endswith('.gz') else open
    with opener(path, 'wb') as fh:
        for name, seq, qual in records:
            fh.write(b'@' + name.encode() + b'\n')
            fh.write(seq + b'\n')
            fh.write(b'+\n')
            fh.write(qual + b'\n')


def parse_vcf(path: str):
    records = {}
    opener = gzip.open if path.endswith('.gz') else open
    with opener(path, 'rt') as fh:
        for line in fh:
            if line.startswith('#'):
                continue
            parts = line.rstrip('\n').split('\t')
            chrom, pos, _, ref, alt, _, _, _, fmt, sample = parts[:10]
            pos = int(pos)
            fmt_fields = fmt.split(':')
            sample_fields = sample.split(':')
            fmt_dict = dict(zip(fmt_fields, sample_fields))
            records[(chrom, pos)] = {
                'ref': ref.encode(),
                'alt': alt,
                'gt': fmt_dict.get('GT', './.'),
                'dp': int(fmt_dict.get('DP', '0')),
                'ad': fmt_dict.get('AD', '0,0'),
            }
    return records


def evaluate(truth, vcf, loci):
    n_truth = len(truth)
    n_called = 0
    correct_gt = 0
    confusion = defaultdict(int)
    alt_tp = alt_fp = alt_fn = 0
    dosage_truth = []
    dosage_call = []

    for (locus_idx, tag_pos), t in truth.items():
        key = (t['contig'], t['genomic_pos'] + 1)  # VCF is 1-based
        t_gt = t['gt']
        t_dosage = {'0/0': 0.0, '0/1': 1.0, '1/1': 2.0}.get(t_gt, float('nan'))
        if key in vcf:
            n_called += 1
            c = vcf[key]
            c_gt = c['gt']
            c_dosage = {'0/0': 0.0, '0/1': 1.0, '1/1': 2.0}.get(c_gt, float('nan'))
            confusion[(t_gt, c_gt)] += 1
            if t_gt == c_gt:
                correct_gt += 1
            # alt allele detection
            if t_gt in ('0/1', '1/1'):
                if c_gt in ('0/1', '1/1'):
                    alt_tp += 1
                else:
                    alt_fn += 1
            else:
                if c_gt in ('0/1', '1/1'):
                    alt_fp += 1
            dosage_truth.append(t_dosage)
            dosage_call.append(c_dosage)
        else:
            confusion[(t_gt, 'NOT_CALLED')] += 1
            if t_gt in ('0/1', '1/1'):
                alt_fn += 1
            dosage_truth.append(t_dosage)
            dosage_call.append(float('nan'))

    concordance = correct_gt / n_truth if n_truth else 0.0
    call_rate = n_called / n_truth if n_truth else 0.0
    alt_recall = alt_tp / (alt_tp + alt_fn) if (alt_tp + alt_fn) else 0.0
    alt_precision = alt_tp / (alt_tp + alt_fp) if (alt_tp + alt_fp) else 0.0
    alt_f1 = 2 * alt_precision * alt_recall / (alt_precision + alt_recall) if (alt_precision + alt_recall) else 0.0

    # dosage correlation ignoring uncalled
    paired = [(t, c) for t, c in zip(dosage_truth, dosage_call) if not (math.isnan(t) or math.isnan(c))]
    if len(paired) >= 2:
        import statistics
        t_mean = statistics.mean(t for t, _ in paired)
        c_mean = statistics.mean(c for _, c in paired)
        num = sum((t - t_mean) * (c - c_mean) for t, c in paired)
        den_t = sum((t - t_mean) ** 2 for t, _ in paired) ** 0.5
        den_c = sum((c - c_mean) ** 2 for _, c in paired) ** 0.5
        dosage_r = num / (den_t * den_c) if den_t and den_c else 0.0
    else:
        dosage_r = 0.0

    return {
        'n_truth': n_truth,
        'n_called': n_called,
        'call_rate': call_rate,
        'concordance': concordance,
        'alt_tp': alt_tp,
        'alt_fp': alt_fp,
        'alt_fn': alt_fn,
        'alt_precision': alt_precision,
        'alt_recall': alt_recall,
        'alt_f1': alt_f1,
        'dosage_r': dosage_r,
        'confusion': dict(confusion),
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--host-db', required=True)
    parser.add_argument('--enzyme', default='BcgI')
    parser.add_argument('--n-loci', type=int, default=1000)
    parser.add_argument('--coverage', type=int, default=20)
    parser.add_argument('--error-rate', type=float, default=0.05)
    parser.add_argument('--error-qual', type=int, default=20,
                        help='Phred quality assigned to simulated errors (default 20 passes -q 20)')
    parser.add_argument('--seed', type=int, default=42)
    parser.add_argument('--min-depth', type=int, default=4)
    parser.add_argument('--min-qual', type=int, default=20)
    parser.add_argument('--max-mismatch', type=int, default=2)
    parser.add_argument('--threads', type=int, default=4)
    parser.add_argument('--outdir', required=True)
    parser.add_argument('--bin', default='f2brad-host')
    args = parser.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    print(f"Loading first {args.n_loci} loci from {args.host_db}")
    loci = load_host_db(args.host_db, args.n_loci)
    print(f"Simulating {args.n_loci} loci at {args.coverage}x with error rate {args.error_rate} (error Q{args.error_qual})")
    truth, r1, r2 = simulate(loci, args.coverage, args.error_rate, args.error_qual, args.seed)

    r1_path = outdir / 'sim_R1.fq.gz'
    r2_path = outdir / 'sim_R2.fq.gz'
    write_fastq(str(r1_path), r1)
    write_fastq(str(r2_path), r2)

    truth_path = outdir / 'truth.tsv'
    with open(truth_path, 'w') as fh:
        fh.write('locus_idx\tcontig\tgenomic_pos\tref\talt\tgt\n')
        for (locus_idx, tag_pos), t in truth.items():
            fh.write(f"{locus_idx}\t{t['contig']}\t{t['genomic_pos']}\t{t['ref'].decode()}\t{t['alt'].decode()}\t{t['gt']}\n")

    gt_out = outdir / 'genotypes'
    gt_out.mkdir(exist_ok=True)
    cmd = [
        args.bin, 'genotype',
        '-d', args.host_db,
        '-1', str(r1_path),
        '-2', str(r2_path),
        '-s', args.enzyme,
        '-o', str(gt_out),
        '-q', str(args.min_qual),
        '--min-depth', str(args.min_depth),
        '--max-mismatch', str(args.max_mismatch),
        '-j', str(args.threads),
    ]
    print('Running:', ' '.join(cmd))
    subprocess.run(cmd, check=True)

    vcf = parse_vcf(str(gt_out / 'genotypes.vcf'))
    metrics = evaluate(truth, vcf, loci)

    print('\n=== Results ===')
    for k, v in metrics.items():
        if k == 'confusion':
            print('confusion matrix (truth_gt -> called_gt):')
            for (tg, cg), cnt in sorted(v.items()):
                print(f"  {tg} -> {cg}: {cnt}")
        else:
            print(f"{k}: {v}")

    with open(outdir / 'metrics.json', 'w') as fh:
        import json
        json_metrics = dict(metrics)
        json_metrics['confusion'] = {f"{tg}->{cg}": cnt for (tg, cg), cnt in metrics['confusion'].items()}
        json.dump(json_metrics, fh, indent=2)


if __name__ == '__main__':
    main()
