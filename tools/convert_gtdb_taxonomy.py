#!/usr/bin/env python3
"""
convert_gtdb_taxonomy.py

功能:
- 将 GTDB 三列表 (accession, gtdb_taxonomy, [ncbi_taxonomy]) 转换为 fast2bRAD-M/2bRAD-M quantify 期望的分类文件:
  abfh_classify_with_speciename.txt.gz

输出格式(制表符分隔):
GCF_or_GCA_ID  Kingdom  Phylum  Class  Order  Family  Genus  Species  Strain

使用:
  python tools/convert_gtdb_taxonomy.py \
    --input /abs/path/test_genomes_taxonomy.tsv \
    --output-dir /abs/path/out_dir
生成:
  /abs/path/out_dir/abfh_classify_with_speciename.txt.gz
"""
import argparse
import gzip
import os
import sys


def extract_gcf_id(filename: str) -> str:
    """
    从文件名式 accession 提取标准 GCF/GCA ID。
    规则与 Rust 版 build_* 保持一致:
    - 去掉 _genomic 及其后缀
    - 按 '_' 切分, 取前两个段, 组合成 GCF_XXXXXXXXX.Y / GCA_XXXXXXXXX.Y
    """
    name = filename.split('/')[-1]
    pos = name.find("_genomic")
    if pos != -1:
        name_clean = name[:pos]
    else:
        name_clean = name

    if name_clean.startswith("GCF_") or name_clean.startswith("GCA_"):
        parts = name_clean.split('_')
        if len(parts) >= 2:
            return f"{parts[0]}_{parts[1]}"
    return name_clean


def parse_gtdb_taxonomy(gtdb_str: str, genome_id: str):
    """
    解析 GTDB taxonomy: d__..;p__..;...;s__..
    返回 8 级: kingdom..strain
    未提供真实株系时, 第 8 级 strain = species + ' ' + genome_id, 使每个基因组成为独立株系。
    必须与 Rust build_qual_db.rs / build_quan_db.rs 的 parse_gtdb_taxonomy 保持一致,
    否则 strain 级 quantify 会按本文件的 Strain 列把同物种基因组重新合并。
    """
    parts = [p.strip() for p in gtdb_str.split(';') if p.strip()]
    if len(parts) < 7:
        raise ValueError(f"GTDB 分类格式错误, 需要至少7级: {gtdb_str}")

    def strip_prefix(x: str) -> str:
        if '__' in x:
            return x.split('__', 1)[1]
        return x

    vals = [strip_prefix(x) for x in parts[:7]]
    species = vals[6]
    strain = f"{species} {genome_id}"
    vals.append(strain)
    return vals  # 8 级


def main():
    ap = argparse.ArgumentParser(description="Convert GTDB 3-col taxonomy TSV to abfh_classify_with_speciename.txt.gz")
    ap.add_argument("--input", "-i", required=True, help="GTDB 三列表 TSV (accession, gtdb_taxonomy, [ncbi_taxonomy])")
    ap.add_argument("--output-dir", "-o", required=True, help="输出目录, 将生成 abfh_classify_with_speciename.txt.gz")
    args = ap.parse_args()

    os.makedirs(args.output_dir, exist_ok=True)
    out_path = os.path.join(args.output_dir, "abfh_classify_with_speciename.txt.gz")

    count = 0
    with open(args.input, "r", encoding="utf-8") as fin, gzip.open(out_path, "wt", encoding="utf-8") as fout:
        for line in fin:
            line = line.strip()
            if not line or line.startswith("#") or line.startswith("accession"):
                continue
            cols = line.split("\t")
            if len(cols) < 2:
                continue
            gcf = extract_gcf_id(cols[0])
            try:
                tax = parse_gtdb_taxonomy(cols[1], gcf)
            except Exception as e:
                print(f"[WARN] 跳过行(分类解析失败): {line}\n  原因: {e}", file=sys.stderr)
                continue
            fout.write(gcf + "\t" + "\t".join(tax) + "\n")
            count += 1

    print(f"[OK] 已写出 {count} 条记录 -> {out_path}")


if __name__ == "__main__":
    main()


