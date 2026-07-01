# 2bRAD-M vs fast2bRAD-M 功能对比

## 功能模块对比

### ✅ 已完全移植的功能

| 2bRAD-M 脚本 | fast2bRAD-M 命令 | 状态 | 说明 |
|-------------|-----------------|------|------|
| `2bRADExtraction.pl` | `extract` | ✅ 完成 | 支持 Type 1-4 全部输入类型，输出与 Perl 版本一致 |
| `CreateQualDatabase_2bRAD.pl` | `build-qual-db` | ✅ 完成 | 构建定性数据库，输出所有标签+unique标记 |
| `CreateQuanDatabase_2bRAD.pl` | `build-quan-db` | ✅ 完成 | 构建定量数据库，只输出unique标签 |
| `CalculateRelativeAbundance_Single2bEnzyme.pl` | `quantify` | ✅ 完成 | 单酶丰度计算，可用于定性和定量分析 |
| `MergeProfilesFromMultipleSamples.pl` | `merge` | ✅ 完成 | 合并多样品丰度表，支持mock和control过滤 |
| `FindGenome_ByQualitative.pl` | `find-genome` | ✅ 完成 | 根据定性结果筛选定量基因组，输出与 Perl 版本一致 |
| `CalculateRelativeAbundance_Combined2bEnzymes.pl` | `quantify` | ⚠️ 部分 | 多酶结果合并逻辑未单独实现，需要手动多次调用quantify后合并 |

### ✅ 已实现的功能（新增）

| 功能 | fast2bRAD-M 命令 | 状态 | 说明 |
|------|-----------------|------|------|
| 一体化主流程 | `pipeline` | ✅ 完成 | `pipeline` 子命令自动串联所有步骤 |
| 机器学习污染分类 | `classify` | ✅ 完成 | 基于 ONNX 模型的污染检测，输出预测标签 |
| 功能预测 | `predict` | ✅ 完成 | 矩阵乘法计算功能丰度（KO/KEGG） |

### ❌ 未实现的功能

| 2bRAD-M 脚本 | 功能描述 | 状态 | 影响 |
|-------------|---------|------|------|
| `2bRADM_Pipline.pl` | 一体化主流程脚本 | ✅ 已替代 | `pipeline` 子命令已完全替代 |

## 详细功能分析

### 1. 数字酶切 (extract) ✅

**2bRAD-M**: `scripts/2bRADExtraction.pl`
- 支持 Type 1-4 全部输入类型
- 支持 16 种酶
- 支持质量控制

**fast2bRAD-M**: `extract` 子命令
- ✅ 完全实现
- ✅ 支持批量并行处理
- ✅ 输出格式与 Perl 版本一致

### 2. 构建定性数据库 (build-qual-db) ✅

**2bRAD-M**: `scripts/CreateQualDatabase_2bRAD.pl`
- 从参考基因组构建分类特异性标签数据库
- 输出所有标签+unique标记

**fast2bRAD-M**: `build-qual-db` 子命令
- ✅ 完全实现
- ✅ 支持多分类层级
- ✅ 输出格式一致

### 3. 构建定量数据库 (build-quan-db) ✅

**2bRAD-M**: `scripts/CreateQuanDatabase_2bRAD.pl`
- 只输出unique标签
- 用于定量分析

**fast2bRAD-M**: `build-quan-db` 子命令
- ✅ 完全实现
- ✅ 支持预酶切文件输入
- ✅ 输出格式一致

### 4. 定性/定量分析 (quantify) ✅

**2bRAD-M**: 
- `CalculateRelativeAbundance_Single2bEnzyme.pl` - 单酶分析
- `CalculateRelativeAbundance_Combined2bEnzymes.pl` - 多酶合并

**fast2bRAD-M**: `quantify` 子命令
- ✅ 单酶分析完全实现
- ⚠️ 多酶合并需要手动多次调用后合并结果

**使用方式对比**：
```bash
# 2bRAD-M: 单酶分析
perl CalculateRelativeAbundance_Single2bEnzyme.pl -l list -d db -t species -s 5 -o out

# 2bRAD-M: 多酶合并
perl CalculateRelativeAbundance_Combined2bEnzymes.pl -l list -s 5,6,7 -io out -m combine

# fast2bRAD-M: 单酶分析（与Perl版本等价）
fast2bRAD-M quantify -l list -d db -t species -s 5 -o out

# fast2bRAD-M: 多酶分析（需要多次调用）
fast2bRAD-M quantify -l list -d db -t species -s 5 -o out
fast2bRAD-M quantify -l list -d db -t species -s 6 -o out
fast2bRAD-M quantify -l list -d db -t species -s 7 -o out
# 然后手动合并结果
```

### 5. 根据定性结果找基因组 (FindGenome) ✅

**2bRAD-M**: `scripts/FindGenome_ByQualitative.pl`

**fast2bRAD-M**: `find-genome` 子命令
- ✅ 完全实现
- ✅ 自动根据 G-score 和 GCF 阈值筛选
- ✅ 输出 `sdb.list` 用于定量建库

```bash
fast2bRAD-M find-genome \
  -l sample_list.tsv \
  -d database_dir \
  -o output_dir \
  --qual-dir qualitative_dir \
  --gscore 5 \
  --gcf 1
```

### 6. 结果合并 (merge) ✅

**2bRAD-M**: `scripts/MergeProfilesFromMultipleSamples.pl`

**fast2bRAD-M**: `merge` 子命令
- ✅ 完全实现
- ✅ 支持 mock 样品过滤
- ✅ 支持阴性对照过滤
- ✅ 输出 `all.xls` 和 `filtered.xls`

### 7. 功能预测 (predict) ✅

**fast2bRAD-M**: `predict` 子命令
- ✅ 矩阵乘法计算功能丰度
- ✅ 输入物种丰度表 + 物种→功能映射矩阵
- ✅ 输出 `.func.xls`

### 8. 机器学习污染分类 (classify) ✅ 新增

**fast2bRAD-M**: `classify` 子命令
- ✅ 基于 ONNX 模型的污染检测
- ✅ 输入 `quantify` 输出文件，输出带 `Prediction` 标签
- ✅ 使用 4 个特征：覆盖率、G-score、深度、理论reads比例

```bash
fast2bRAD-M classify \
  -i sample.BcgI.xls \
  -m contamination_model.onnx \
  -o sample.BcgI.classified.xls
```

### 9. 一体化主流程 (Pipeline) ✅ 新增

**fast2bRAD-M**: `pipeline` 子命令
- ✅ 自动串联所有步骤
- ✅ 支持 `full` / `db-only` / `sample-only` 三种模式
- ✅ 支持断点续跑（`.done` 标记文件）
- ✅ 支持功能预测（`--ko-mapping`）

```bash
# 完整流程
fast2bRAD-M pipeline \
  --mode full \
  --samples samples.tsv \
  --genome-list genome_list.tsv \
  --taxonomy taxonomy.tsv \
  --site BcgI --level species \
  --outdir results/ --prefix run1 --threads 16

# 仅建库
fast2bRAD-M pipeline --mode db-only \
  --genome-list genome_list.tsv \
  --taxonomy taxonomy.tsv \
  --site BcgI --level species \
  --outdir db/ --threads 16

# 仅分析（使用已有数据库）
fast2bRAD-M pipeline --mode sample-only \
  --samples samples.tsv \
  --database db/ \
  --site BcgI --level species \
  --outdir results/ --prefix run1 --threads 16
```

## 缺失功能影响评估

### 高优先级缺失
无 — 核心功能已完全覆盖。

### 中优先级

1. **CalculateRelativeAbundance_Combined2bEnzymes.pl** ⚠️
   - **影响**：多酶分析需要手动合并结果
   - **解决方案**：可以多次调用 `quantify` 后手动合并
   - **建议**：在 `quantify` 中增加多酶支持，或添加 `combine` 子命令

## 功能完整性总结

### 核心功能完整性：~98%

- ✅ 数字酶切：100%
- ✅ 数据库构建：100%
- ✅ 单酶分析：100%
- ✅ FindGenome：100%
- ✅ 结果合并：100%
- ✅ 功能预测：100%
- ✅ 机器学习分类：100%
- ✅ 一体化流程：100%
- ⚠️ 多酶合并：50%（需要手动）

### 建议实现顺序

1. ~~实现 `find-genome` 子命令~~ ✅ 已完成
2. ~~实现 `pipeline` 子命令~~ ✅ 已完成
3. ~~实现 `classify` 子命令（ML 污染分类）~~ ✅ 已完成
4. **中优先级**：增强 `quantify` 支持多酶，或添加 `combine` 子命令
5. **低优先级**：性能 benchmark 对比 Perl 版本

## 使用建议

### 单酶分析（完整流程）

fast2bRAD-M 已完全支持单酶分析的完整流程：

```bash
# 方式1：一键 pipeline
fast2bRAD-M pipeline \
  --mode full \
  --samples samples.tsv \
  --genome-list genome_list.tsv \
  --taxonomy taxonomy.tsv \
  --site BcgI --level species \
  --outdir results/ --prefix run1 --threads 16

# 方式2：分步执行（更灵活）
fast2bRAD-M extract --batch samples.tsv -t 2 -s BcgI --od 01_extract
fast2bRAD-M build-qual-db -l genome_list.tsv -s BcgI -t species -o 02_db_qual
fast2bRAD-M quantify -l samples.tsv -d 02_db_qual -t species -s BcgI -o qualitative
fast2bRAD-M find-genome -l samples.tsv -d 02_db_qual -o sdb --qual-dir qualitative
fast2bRAD-M build-quan-db -l sdb.list -s BcgI -t species -o 02_db_quan
fast2bRAD-M quantify -l samples.tsv -d 02_db_quan -t species -s BcgI -o 04_quantify
fast2bRAD-M merge -l merge_list.tsv -o 05_merge -p Abundance_Stat

# 可选：功能预测
fast2bRAD-M predict -a 05_merge/Abundance_Stat.all.xls -m ko_mapping.tsv -o 05_merge

# 可选：ML 污染分类
fast2bRAD-M classify -i 04_quantify/sample1/sample1.BcgI.xls -m model.onnx -o sample1.classified.xls
```

### 多酶分析

对于多酶分析，当前需要手动执行多次后合并。建议使用脚本组合：

```bash
for enzyme in BcgI CjeI BsaXI; do
  fast2bRAD-M pipeline --mode sample-only \
    --samples samples.tsv \
    --database db/ \
    --site $enzyme --level species \
    --outdir results_${enzyme}/ --prefix ${enzyme} \
    --threads 16
done
# 然后手动合并多个酶的结果
```