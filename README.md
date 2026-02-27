# fast2bRAD-M

高性能 Rust 版 2bRAD-M 微生物组分析工具

## 特性

- ⚡ **极致性能**：Rust 实现 + Rayon 多核并行，批量处理 15 个基因组 < 0.12 秒
- 🚀 **批量并行**：自动利用所有 CPU 核心，显著加速大规模数据处理
- 🧬 **完整酶支持**：支持全部 16 种 Type IIB 限制酶
- 📊 **全部输入类型**：参考基因组、Shotgun、单标签
- 🔬 **质量控制**：N 比例、质量分数、质量百分比检查
- 💾 **灵活输出**：FASTA/FASTQ 格式、gzip 压缩

## 安装 

### 方式一：使用 Conda 环境（推荐）

使用提供的 conda 环境配置文件快速搭建环境：

```bash
conda env create -f fast2bRAD-M/fast2brad_m_conda.yaml -n fast2brad
conda activate fast2brad
cd fast2bRAD-M
cargo build --release
```

环境包含：
- Rust 编译工具链
- zlib、openssl（依赖库）

### 方式二：直接编译

```bash
cd fast2bRAD-M
cargo build --release
```

**注意**：如果使用 Type 2 双端数据或 Type 4 功能，需要单独安装 PEAR：
```bash
conda install -c bioconda pear
```

## 使用

### pipeline - 一键主流程（对齐 Perl 管线）

将 Perl 的 `2bRADM_Pipline.pl` 流程一键化，串联：
extract → （可选）build-qual-db → （可选）build-quan-db → quantify → merge  
产物目录结构：`01_extract/ 02_db_qual/ 04_quantify/ 05_merge/ qualitative/ quantitative_sdb/`

```bash
# 全流程（含数据库构建）：
fast2bRAD-M pipeline \
  --mode full \
  --samples /abs/samples_list.tsv \
  --taxonomy /abs/abfh_classify_with_speciename.txt \
  --pre-digested-dir /abs/pre_digested_output \
  --site BcgI \
  --level species \
  --outdir /abs/outdir \
  --prefix run \
  --gscore 5 \
  --threads 4 \
  --pc 8\
  --min-qual 15 \
  --resume yes

# 数据库构建：
fast2bRAD-M extract \
  --genome-list /abs/pre_digested_file_list.tsv \
  -t 1 \
  -s BcgI \
  --od /abs/pre_digested_output \
  --op db \
  --threads 5 \
  --qc yes \
  -n 0.08 \
  -q 30 \
  -p 80

fast2bRAD-M pipeline \
  --mode db-only \
  --genome-list /abs/pre_digested_file_list.tsv \
  --taxonomy /abs/abfh_classify_with_speciename.txt \
  --site BcgI \
  --level species \
  --outdir /abs/outdir \
  --prefix run \
  --gscore 5 \
  --threads 4 \
  --pc 8 \
  --min-qual 15 \
  --resume yes

# 使用已有数据库（推荐：与 Perl 一致的 classify 与 *.fa.gz）
fast2bRAD-M pipeline \
  --mode sample-only \
  --samples /abs/samples_list.tsv \
  --database /abs/db_ready \
  --site BcgI \
  --level species \
  --outdir /abs/outdir \
  --prefix run \
  --gscore 5 \
  --threads 4 \
  --pc 8 \
  --min-qual 15 \
  --resume yes
```

参数与默认（对齐 Perl 取舍）：
- `--mode`：full|db-only|sample-only（默认 full）
- `--gscore`：默认 5（与 Perl 常用阈值一致）
- `--resume`：默认 yes（存在产物则跳过）
- `--threads`：设置 `RAYON_NUM_THREADS`；不设则自动
- `--mock`、`--control`：合并阶段过滤（与 Perl 行为一致）
- `--samples`：TSV：`sample<TAB>path1[<TAB>path2]`（原始 FASTQ/FASTA 路径，非 .iibsp）
- 数据库目录需包含：`BcgI.species.fa.gz` 和 `abfh_classify_with_speciename.txt.gz`

Bash 包装脚本（与上完全等价）：
```bash
bash fast2bRAD-M/scripts/run_pipeline.sh --help
```

**样品列表格式** (`sample_list.tsv`)：
- 第1列：样品名（输出文件前缀）
- 第2列：输入文件路径
- 第3列：输入文件路径2（可选，用于 Type 2/4 双端数据）
- 以 `#` 开头的行为注释

示例：
```tsv
# Type 1: 参考基因组
ecoli	/path/to/ecoli.fna.gz
lplantarum	/path/to/lplantarum.fna.gz
sagalactiae	/path/to/sagalactiae.fna.gz
```

```tsv
# Type 2: Shotgun 双端测序
sample1	/path/to/sample1_R1.fq.gz	/path/to/sample1_R2.fq.gz
sample2	/path/to/sample2_R1.fq.gz	/path/to/sample2_R2.fq.gz
```

### 输入类型说明

1. **Type 1**: 参考基因组 FASTA - 滑动窗口全匹配，用于构建数据库
2. **Type 2**: Shotgun 测序（SE/PE）- 序列内匹配所有标签位点（去重）
3. **Type 3**: 2bRAD 单标签（SE）- 只取第一个匹配的标签

### 支持的酶（1-16）

1. CspCI   2. AloI     3. BsaXI    4. BaeI  
5. BcgI（默认） 6. CjeI  7. PpiI     8. PsrI  
9. BplI    10. FalI    11. Bsp24I  12. HaeIV  
13. CjePI  14. Hin4I   15. AlfI    16. BslFI

### 参数说明

- `--batch`：批量处理样品列表（TSV 格式）
- `-i, --input`：输入文件（支持 .gz）
- `-t, --type`：输入类型（1-3）
- `-s, --site`：酶编号（1-16）或名称
- `--od`：输出目录
- `--op`：输出前缀
- `--gz`：是否压缩（yes/no，默认 yes）
- `--qc`：是否质控（yes/no，默认 yes）
- `-n, --max-n`：最大 N 比例（默认 0.08）
- `-q, --min-quality`：最低质量分数（默认 30）
- `-p, --min-quality-percent`：最低质量百分比（默认 80）
- `-b, --quality-base`：质量编码（默认 33）
- `--fm`：输出格式（fa/fq，默认 fa）

## 输出

### 标签文件

- Type 1-3 文件名：`{sample}.{enzyme}.{format}.gz`
- 扩展名：`.iibsp`（IIB = Type IIB 限制酶）
- 格式：FASTA 或 FASTQ

### 统计文件

- Type 1-3：`{sample}.{enzyme}.stat.tsv`
- 内容：输入序列数、标签数、百分比等


### build-qual-db - 构建定性数据库

从参考基因组构建分类特异性 2bRAD 标签数据库：

```bash
# 构建 species 级别数据库
fast2bRAD-M build-qual-db \
  -l genome_list.tsv \
  -s 5 \
  -t species \
  -o database_dir

# 构建多个级别数据库
fast2bRAD-M build-qual-db \
  -l genome_list.tsv \
  -s BcgI \
  -t genus,species \
  -o database_dir

# 构建所有级别数据库
fast2bRAD-M build-qual-db \
  -l genome_list.tsv \
  -s 5 \
  -t all \
  -o database_dir \
  -r yes
```

### 基因组列表格式

TSV 文件，每行包含：
```
GCFid    Kingdom    Phylum    Class    Order    Family    Genus    Species    Strain    genome_path
```

示例：
```
GCF_000007445.1    Bacteria    Proteobacteria    Gammaproteobacteria    Enterobacterales    Enterobacteriaceae    Escherichia    Escherichia_coli    str._K-12    /path/to/genome.fna.gz
```

## 已完成功能

✅ **extract**（数字酶切）- 支持 Type 1-3 全部输入类型  
✅ **build-qual-db**（定性数据库）- 输出所有标签+unique标记  
✅ **build-quan-db**（定量数据库）- 只输出unique标签  
✅ **quantify**（丰度计算）- 计算样品中微生物相对丰度，输出 GCF_detected.xls  
✅ **find-genome**（筛选基因组）- 根据定性结果筛选定量分析所需的基因组  
✅ **merge**（结果合并）- 合并多样品丰度表  

## 使用示例

### 完整分析流程

```bash
# 1. 数字酶切
fast2bRAD-M extract --batch sample_list.tsv -t 2 -s 5 --od enzyme_result

# 2. 定性分析
fast2bRAD-M quantify -l enzyme_list.tsv -d qual_db -t species -s 5 -o qualitative

# 3. 根据定性结果筛选基因组（新增功能）
fast2bRAD-M find-genome \
  -l sample_list.tsv \
  -d database_dir \
  -o quantitative_sdb \
  --qual-dir qualitative \
  --gscore 5 \
  --gcf 1

# 4. 构建样品特异性定量数据库
fast2bRAD-M build-quan-db \
  -l quantitative_sdb/sample1/sdb.list \
  -s 5 \
  -t species \
  -o quantitative_sdb/sample1/database

# 5. 定量分析
fast2bRAD-M quantify -l sample_list.tsv -d quantitative_sdb/sample1/database -t species -s 5 -o quantitative

# 6. 结果合并
fast2bRAD-M merge -l abundance_list.tsv -o quantitative -p Abundance_Stat
```

### find-genome - 根据定性结果筛选基因组

根据定性分析结果，筛选出用于定量分析的候选基因组：

```bash
fast2bRAD-M find-genome \
  -l sample_list.tsv \
  -d database_dir \
  -o output_dir \
  --qual-dir qualitative_dir \
  --gscore 5 \
  --gcf 1
```

**参数说明**：
- `-l, --list`: 样品列表文件（TSV格式）
- `-d, --database`: 数据库目录（需包含 `abfh_classify_with_speciename.txt.gz`）
- `-o, --output`: 输出目录
- `--qual-dir`: 定性分析结果目录
- `--gscore`: G-score 阈值（默认 5，表示 >5）
- `--gcf`: GCF 标签数阈值（默认 1，表示 >1）

**输出**：
- `$output_dir/$sample/sdb.list` - 每个样品的候选基因组列表


## pipeline：双端与 PEAR 透传

流水线使用已有数据库时的一键示例（双端+PEAR）：

```bash
fast2bRAD-M pipeline \
  --mode sample-only \
  -l /abs/samples.tsv \
  -s BcgI -t species \
  --outdir /abs/runs/run_pe \
  --prefix run_pe \
  -d /abs/db_ready \
  --pe pear --pc 8 \
  --resume yes
```

说明：
- 若样品行提供两列路径（R1、R2），pipeline 会透传 `--pe/--pc` 到 extract，先调用 PEAR 拼接，再继续提取。
- 生成的中间合并文件为 `<prefix>.<enzyme>.pear.fastq`，最终样品标签文件为 `<prefix>.<enzyme>.iibsp[.gz]`。

## 许可证

继承原 2bRAD-M 项目许可证

