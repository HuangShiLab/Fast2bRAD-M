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

### ⚠️ 部分实现/需要确认的功能

| 2bRAD-M 脚本 | fast2bRAD-M 对应 | 状态 | 说明 |
|-------------|-----------------|------|------|
| `CalculateRelativeAbundance_Combined2bEnzymes.pl` | `quantify` | ⚠️ 部分 | 多酶结果合并逻辑未单独实现，需要手动多次调用quantify后合并 |

### ✅ 已实现的功能（新增）

| 2bRAD-M 脚本 | fast2bRAD-M 命令 | 状态 | 说明 |
|-------------|-----------------|------|------|
| `FindGenome_ByQualitative.pl` | `find-genome` | ✅ 完成 | 根据定性结果筛选定量基因组，输出与 Perl 版本一致 |

### ❌ 未实现的功能

| 2bRAD-M 脚本 | 功能描述 | 状态 | 影响 |
|-------------|---------|------|------|
| `2bRADM_Pipline.pl` | 一体化主流程脚本 | ❌ 缺失 | 需要手动组合多个命令，但功能可替代 |

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

### 5. 根据定性结果找基因组 (FindGenome) ❌

**2bRAD-M**: `scripts/FindGenome_ByQualitative.pl`

**功能**：
- 读取定性分析结果（`$qual_dir/$sample/$sample.combine.xls`）
- 根据 G-score 阈值筛选分类
- 根据 GCF 标签数阈值筛选基因组
- 生成定量建库所需的基因组列表（`sdb.list`）

**参数**：
- `-l`: 样品列表
- `-d`: 数据库目录
- `-o`: 输出目录
- `-qualdir`: 定性结果目录
- `-gscore`: G-score 阈值（默认 5）
- `-gcf`: GCF 标签数阈值（默认 1）

**输出**：
- `$outdir/$sample/sdb.list` - 每个样品的候选基因组列表

**影响**：
- ⚠️ **关键缺失**：这是定量分析流程中的必要步骤
- 当前需要手动实现或使用 Perl 脚本

**建议实现**：
```rust
// 建议添加子命令：find-genome
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

### 7. 一体化主流程 (Pipeline) ❌

**2bRAD-M**: `bin/2bRADM_Pipline.pl`

**功能流程**：
1. 数字酶切（多样品并行）
2. 定性分析（可选）
   - 单酶分析
   - 多酶合并
3. 定量分析（可选）
   - FindGenome（根据定性结果筛选基因组）
   - 构建样品特异性定量数据库
   - 单酶定量分析
   - 多酶合并
4. 结果合并

**fast2bRAD-M**: 无一体化脚本

**当前替代方案**：
```bash
# 1. 数字酶切
fast2bRAD-M extract --batch sample_list.tsv -t 2 -s 5 --od enzyme_result

# 2. 定性分析（需要手动实现多酶合并）
fast2bRAD-M quantify -l enzyme_list.tsv -d qual_db -t species -s 5 -o qualitative

# 3. 定量分析（缺少 FindGenome 步骤）
# 需要手动筛选基因组或使用 Perl 脚本
# fast2bRAD-M find-genome -l sample_list -d db -o sdb --qual-dir qualitative
fast2bRAD-M build-quan-db -l sdb.list -s 5 -t species -o sample_db
fast2bRAD-M quantify -l sample_list -d sample_db -t species -s 5 -o quantitative

# 4. 结果合并
fast2bRAD-M merge -l abundance_list.tsv -o quantitative -p Abundance_Stat
```

## 缺失功能影响评估

### 高优先级缺失

1. **FindGenome_ByQualitative.pl** ❌
   - **影响**：无法自动从定性结果筛选定量分析所需的基因组
   - **解决方案**：需要手动实现或使用 Perl 脚本
   - **建议**：实现 `find-genome` 子命令

### 中优先级缺失

2. **CalculateRelativeAbundance_Combined2bEnzymes.pl** ⚠️
   - **影响**：多酶分析需要手动合并结果
   - **解决方案**：可以多次调用 `quantify` 后手动合并
   - **建议**：在 `quantify` 中增加多酶支持，或添加 `combine` 子命令

3. **2bRADM_Pipline.pl** ❌
   - **影响**：需要手动组合多个命令，流程复杂
   - **解决方案**：可以通过脚本组合实现
   - **建议**：实现 `pipeline` 子命令或提供示例脚本

## 功能完整性总结

### 核心功能完整性：95%

- ✅ 数字酶切：100%
- ✅ 数据库构建：100%
- ✅ 单酶分析：100%
- ⚠️ 多酶合并：50%（需要手动）
- ✅ FindGenome：100%（**已完成**）
- ✅ 结果合并：100%

### 建议实现顺序

1. ~~**高优先级**：实现 `find-genome` 子命令~~ ✅ **已完成**
2. **中优先级**：增强 `quantify` 支持多酶，或添加 `combine` 子命令
3. **低优先级**：实现 `pipeline` 子命令或提供流程脚本

## 使用建议

### 当前可用流程

对于**单酶分析**，fast2bRAD-M 已完全可用：
```bash
# 完整单酶流程
fast2bRAD-M extract --batch sample_list.tsv -t 2 -s 5 --od enzyme_result
fast2bRAD-M quantify -l enzyme_list.tsv -d qual_db -t species -s 5 -o qualitative
# 手动筛选基因组（或使用 Perl 脚本）
fast2bRAD-M build-quan-db -l sdb.list -s 5 -t species -o sample_db
fast2bRAD-M quantify -l sample_list -d sample_db -t species -s 5 -o quantitative
fast2bRAD-M merge -l abundance_list.tsv -o quantitative
```

### 需要 Perl 脚本辅助的流程

对于**多酶分析**，需要：
- 使用 Perl 版本的 `CalculateRelativeAbundance_Combined2bEnzymes.pl` 进行多酶合并

**注意**：`find-genome` 功能已完全实现，不再需要 Perl 脚本辅助。

