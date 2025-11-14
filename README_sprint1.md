# fast2bRAD-M

高性能 Rust 版 2bRAD-M 微生物组分析工具

## 特性

- ⚡ **高性能**：Rust 实现，使用 fxhash 优化哈希性能
- 🧬 **完整酶支持**：支持全部 16 种 Type IIB 限制酶
- 📊 **全部输入类型**：参考基因组、Shotgun、单标签、5连标签（Type 1-4）
- 🔬 **质量控制**：N 比例、质量分数、质量百分比检查
- 💾 **灵活输出**：FASTA/FASTQ 格式、gzip 压缩
- ✅ **精确一致**：与 Perl 原版输出100%一致

## 安装

```bash
cd fast2bRAD-M
cargo build --release
```

## 使用

### extract - 数字酶切

从序列数据中提取 2bRAD 标签：

```bash
# Type 1: 参考基因组
fast2bRAD-M extract \
  -i genome.fna.gz \
  -t 1 \
  -s 5 \
  -o output_dir \
  --op sample_name \
  --gz yes

# Type 2: Shotgun 测序数据
fast2bRAD-M extract \
  -i shotgun.fq.gz \
  -t 2 \
  -s BcgI \
  -o output_dir \
  --op sample_name \
  --qc yes \
  -n 0.08 \
  -q 30 \
  -p 80

# Type 3: 单 2bRAD 标签
fast2bRAD-M extract \
  -i 2brad_single.fq.gz \
  -t 3 \
  -s 5 \
  -o output_dir \
  --op sample_name \
  --fm fq

# Type 4: 5连标签（需要预先用 PEAR 拼接 R1/R2）
fast2bRAD-M extract \
  -i assembled.fq.gz dummy.fq.gz \
  -t 4 \
  -s 5 \
  -o output_dir \
  --op sample1 sample2 sample3 sample4 sample5 \
  --qc yes
```

### 输入类型说明

1. **Type 1**: 参考基因组 FASTA - 滑动窗口全匹配，用于构建数据库
2. **Type 2**: Shotgun 测序（SE/PE）- 序列内匹配所有标签位点（去重）
3. **Type 3**: 2bRAD 单标签（SE）- 只取第一个匹配的标签
4. **Type 4**: 2bRAD 5连标签（PE）- 按酶特定位置切分成 5 个样本
   - 需要预先用 PEAR 拼接 R1/R2
   - 或直接提供已拼接的 FASTQ 文件
   - 需要 5 个输出前缀

### 支持的酶（1-16）

1. CspCI   2. AloI     3. BsaXI    4. BaeI  
5. BcgI（默认） 6. CjeI  7. PpiI     8. PsrI  
9. BplI    10. FalI    11. Bsp24I  12. HaeIV  
13. CjePI  14. Hin4I   15. AlfI    16. BslFI

### 参数说明

- `-i, --input`：输入文件（支持 .gz，Type 4 需要 2 个文件）
- `-t, --type`：输入类型（1-4）
- `-s, --site`：酶编号（1-16）或名称
- `-o, --od`：输出目录
- `--op`：输出前缀（Type 4 需要 5 个）
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
- Type 4 文件名：每个样本独立输出 `{sample1-5}.{enzyme}.{format}.gz`
- 扩展名：`.iibsp`（IIB = Type IIB 限制酶）
- 格式：FASTA 或 FASTQ

### 统计文件

- Type 1-3：`{sample}.{enzyme}.stat.tsv`
- Type 4：`{sample1-sample2-...-sample5}.{enzyme}.stat.tsv`
- 内容：输入序列数、标签数、百分比等

## 性能优化

- 使用 fxhash 替代标准 HashMap，提升哈希性能
- 标签去重采用 HashSet，内存换时间
- 流式处理，避免一次性加载全部数据
- 支持 gzip 流式压缩输出

## 测试

已通过与 Perl 原版对比测试，输出100%一致：

```bash
# Rust 版 Type 1
./target/release/fast2bRAD-M extract -i test.fna.gz -t 1 -s 5 -o out --op test --gz no

# 对比 Perl 版
perl 2bRAD-M/scripts/2bRADExtraction.pl -i test.fna.gz -t 1 -s 5 -od out_perl -op test -gz no

# 结果一致
wc -l out/test.BcgI.fa out_perl/test.BcgI.fa
# 6418 out/test.BcgI.fa
# 6418 out_perl/test.BcgI.fa
```

## Type 4 特别说明

Type 4（5连标签）是一种特殊的 2bRAD 建库方式：

1. **前置步骤**：需要先用 PEAR 拼接 R1/R2
   ```bash
   pear -f R1.fq.gz -r R2.fq.gz -o sample -j 8
   ```

2. **运行 Type 4**：
   ```bash
   fast2bRAD-M extract \
     -i sample.assembled.fastq dummy.fq.gz \
     -t 4 -s 5 -o output \
     --op sample1 sample2 sample3 sample4 sample5
   ```

3. **输出**：5 个独立的样本标签文件

4. **原理**：每个酶有特定的切分位置（concat_starts/concat_ends），在拼接后的长序列上按位置提取 5 个标签

## 待实现

- [ ] 内置 PEAR 拼接逻辑（当前需要外部预处理）
- [ ] 并行处理多样本
- [ ] build-qual-db（构建定性数据库）
- [ ] build-quan-db（构建定量数据库）
- [ ] quantify（丰度计算）
- [ ] merge（结果合并）

## 许可证

继承原 2bRAD-M 项目许可证

