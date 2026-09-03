#!/usr/bin/env python3
"""Real-reference-context genotype benchmark.

Extracts genomic fragments around BcgI host tags from the T2T-CHM13 reference,
introduces known SNPs inside the tags, simulates paired-end reads, and compares
the called VCF to the truth.
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


def load_host_db(path: str):
    loci = []
    opener = gzip.open if path.endswith('.gz') else open
    with opener(path, 'rt') as fh:
        reader = csv.DictReader(fh, delimiter='\t')
        for row in reader:
            loci.append({
                'contig': row['contig'],
                'pos': int(row['pos']),
                'strand': row['strand'],
                'seq': row['seq'].encode(),
                'canonical': row['canonical'].encode(),
            })
    return loci


def load_fai(path: str):
    lengths = {}
    with open(path) as fh:
        for line in fh:
            parts = line.split('\t')
            lengths[parts[0]] = int(parts[1])
    return lengths


def add_errors(seq: bytes, error_rate: float, error_qual: int) -> tuple[bytes, bytes]:
    seq = bytearray(seq)
    qual = bytearray()
    for i, base in enumerate(seq):
        if random.random() < error_rate:
            seq[i] = mutate_base(bytes([base]))[0]
            qual.append(error_qual + 33)
        else:
            qual.append(40 + 33)
    return bytes(seq), bytes(qual)


def generate(loci, ref_fa, flank, coverage, error_rate, error_qual, read_len, seed):
    random.seed(seed)
    truth = {}
    r1_records = []
    r2_records = []
    frag_id = 0

    for locus_idx, locus in enumerate(loci):
        tag_len = len(locus['seq'])
        alt_pos = tag_len // 2
        ref_base = bytes([locus['seq'][alt_pos]])
        alt_base = mutate_base(ref_base)
        gt = random.choice(['0/0', '0/1', '1/1'])

        # Genomic coordinate of the SNP on the plus strand.
        if locus['strand'] == '-':
            snp_pos = locus['pos'] + (tag_len - 1 - alt_pos)
        else:
            snp_pos = locus['pos'] + alt_pos

        truth[(locus['contig'], snp_pos + 1)] = {  # 1-based for VCF comparison
            'locus_idx': locus_idx,
            'ref': ref_base,
            'alt': alt_base,
            'gt': gt,
        }

        region_start = max(0, locus['pos'] - flank)
        region_end = locus['pos'] + tag_len + flank
        tag_offset = locus['pos'] - region_start

        for _ in range(coverage):
            is_alt = False
            if gt == '1/1':
                is_alt = True
            elif gt == '0/1':
                is_alt = random.random() < 0.5

            frag = bytearray(locus['fragment'])
            if is_alt:
                frag[tag_offset + alt_pos] = alt_base[0]
            frag = bytes(frag)

            r1_seq, r1_qual = add_errors(frag[:read_len], error_rate, error_qual)
            r2_seq, r2_qual = add_errors(reverse_complement(frag[-read_len:]), error_rate, error_qual)

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
            }
    return records


def evaluate(truth, vcf):
    n = len(truth)
    called = 0
    correct = 0
    alt_tp = alt_fp = alt_fn = 0
    for key, t in truth.items():
        if key in vcf:
            called += 1
            c = vcf[key]
            if t['gt'] == c['gt']:
                correct += 1
            if t['gt'] in ('0/1', '1/1'):
                alt_tp += 1 if c['gt'] in ('0/1', '1/1') else 0
                alt_fn += 1 if c['gt'] not in ('0/1', '1/1') else 0
            else:
                alt_fp += 1 if c['gt'] in ('0/1', '1/1') else 0
        else:
            if t['gt'] in ('0/1', '1/1'):
                alt_fn += 1
    return {
        'n_truth': n,
        'n_called': called,
        'call_rate': called / n,
        'concordance': correct / n,
        'alt_recall': alt_tp / (alt_tp + alt_fn) if (alt_tp + alt_fn) else 0.0,
        'alt_precision': alt_tp / (alt_tp + alt_fp) if (alt_tp + alt_fp) else 0.0,
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--host-db', required=True)
    parser.add_argument('--ref-fa', required=True)
    parser.add_argument('--n-loci', type=int, default=2000)
    parser.add_argument('--flank', type=int, default=150)
    parser.add_argument('--coverage', type=int, default=20)
    parser.add_argument('--error-rate', type=float, default=0.01)
    parser.add_argument('--error-qual', type=int, default=20)
    parser.add_argument('--read-len', type=int, default=150)
    parser.add_argument('--seed', type=int, default=42)
    parser.add_argument('--min-depth', type=int, default=4)
    parser.add_argument('--threads', type=int, default=4)
    parser.add_argument('--outdir', required=True)
    parser.add_argument('--bin', default='f2brad-host')
    args = parser.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    print(f"Loading host DB: {args.host_db}")
    all_loci = load_host_db(args.host_db)
    print(f"Loaded {len(all_loci)} loci")

    fai = load_fai(args.ref_fa + '.fai')

    eligible = []
    for locus in all_loci:
        tag_len = len(locus['seq'])
        length = fai.get(locus['contig'], 0)
        if locus['pos'] >= args.flank and locus['pos'] + tag_len + args.flank <= length:
            eligible.append(locus)

    print(f"Eligible loci with {args.flank} bp flanks: {len(eligible)}")
    random.seed(args.seed)
    selected = eligible[:args.n_loci] if len(eligible) >= args.n_loci else eligible
    print(f"Using first {len(selected)} eligible loci")

    # Extract oriented fragments with bedtools.
    bed_path = outdir / 'regions.bed'
    fa_path = outdir / 'regions.fa'
    with open(bed_path, 'w') as fh:
        for i, locus in enumerate(selected):
            tag_len = len(locus['seq'])
            start = max(0, locus['pos'] - args.flank)
            end = locus['pos'] + tag_len + args.flank
            fh.write(f"{locus['contig']}\t{start}\t{end}\tloc{i}\t0\t{locus['strand']}\n")

    # Extract the plus-strand fragment. The genotyper canonicalizes reads
    # internally, so we do not reverse-complement for minus-strand loci here.
    cmd = ['bedtools', 'getfasta', '-fi', args.ref_fa, '-bed', str(bed_path),
           '-nameOnly', '-fo', str(fa_path)]
    print('Running:', ' '.join(cmd))
    subprocess.run(cmd, check=True)

    # Parse extracted fragments and attach to loci.
    fragments = {}
    current_name = None
    current_seq = []
    with open(fa_path) as fh:
        for line in fh:
            line = line.strip()
            if line.startswith('>'):
                if current_name is not None:
                    fragments[current_name] = ''.join(current_seq).encode()
                current_name = line[1:].split('(')[0]
                current_seq = []
            else:
                current_seq.append(line)
        if current_name is not None:
            fragments[current_name] = ''.join(current_seq).encode()


    for i, locus in enumerate(selected):
        name = f"loc{i}"
        if name not in fragments:
            raise RuntimeError(f"Missing fragment for {name}")
        locus['fragment'] = fragments[name]

    print(f"Simulating {len(selected)} loci at {args.coverage}x")
    truth, r1, r2 = generate(selected, args.ref_fa, args.flank, args.coverage,
                              args.error_rate, args.error_qual, args.read_len, args.seed)

    r1_path = outdir / 'sim_R1.fq.gz'
    r2_path = outdir / 'sim_R2.fq.gz'
    write_fastq(str(r1_path), r1)
    write_fastq(str(r2_path), r2)

    truth_path = outdir / 'truth.tsv'
    with open(truth_path, 'w') as fh:
        fh.write('contig\tpos\tref\talt\tgt\n')
        for (chrom, pos), t in truth.items():
            fh.write(f"{chrom}\t{pos}\t{t['ref'].decode()}\t{t['alt'].decode()}\t{t['gt']}\n")

    gt_out = outdir / 'genotypes'
    gt_out.mkdir(exist_ok=True)
    cmd = [
        args.bin, 'genotype',
        '-d', args.host_db,
        '-1', str(r1_path),
        '-2', str(r2_path),
        '-s', 'BcgI',
        '-o', str(gt_out),
        '--min-depth', str(args.min_depth),
        '-j', str(args.threads),
    ]
    print('Running:', ' '.join(cmd))
    subprocess.run(cmd, check=True)

    vcf = parse_vcf(str(gt_out / 'genotypes.vcf'))
    metrics = evaluate(truth, vcf)

    print('\n=== Results ===')
    for k, v in metrics.items():
        print(f"{k}: {v}")


if __name__ == '__main__':
    main()
