# 实验总结：BRCA1 BRCT 结构域突变效应预测

本文档总结了将项目流程扩展到 BRCA1 C 端 BRCT repeats 区域后的实验结果。由于 BRCA1 全长约 1863 aa，超过当前 ESM-2/SaProt 全长输入的安全长度，本实验不直接对 BRCA1 全长评分，而是采用 domain-local scoring 分析 BRCT 区域，残基位置为 1646-1859。

## 1. 项目目标

本实验的目标是检验在长蛋白的功能结构域局部场景中，结构感知（structure-aware）的蛋白质语言模型和 pLDDT-aware 结构 token 掩蔽是否仍具有价值。

当前核心问题是：

- 在 BRCA1 BRCT 区域的零样本致病性排序中，`SaProt (Seq + 3Di)` 是否优于仅使用序列的 `ESM-2`？
- 基于 pLDDT 对 BRCT 区域中的中低置信度 `3Di` token 进行掩蔽，是否能提升 `SaProt` 表现？
- 对长蛋白进行 domain-local scoring 是否可作为全长输入受限时的合理扩展方案？

## 2. 数据集与预处理

### 目标对象

- 蛋白质：人类 BRCA1
- UniProt ID：`P38398`
- 分析区域：BRCT repeats，残基范围 `1646-1859`
- 结构来源：AlphaFold 模型 `AF-P38398-F1-model_v6.pdb`
- 结构分词（Tokenization）：`Foldseek structureto3didescriptor`
- 评分方式：domain-local scoring，将全长坐标转换为 BRCT 局部坐标后输入模型

### ClinVar 数据清洗规则

数据集基于 ClinVar 的 `variant_summary.txt.gz` 构建，并使用与 TP53/MSH2 实验一致的过滤条件：

- 仅保留 `GeneSymbol == BRCA1` 的记录。
- 仅保留残基位置位于 `1646-1859` 的 BRCT 区域变异。
- 仅解析单氨基酸替换（single amino-acid substitutions）。
- 移除 uncertain、conflicting、not-provided、association-only、drug-response、risk-factor、protective 以及其他非二分类的临床标签。
- 将 `Pathogenic/Likely pathogenic`（致病/疑似致病）合并为标签 `1`。
- 将 `Benign/Likely benign`（良性/疑似良性）合并为标签 `0`。
- 要求 ClinVar 中的野生型残基与 UniProt 参考序列一致。
- 要求野生型残基与从 AlphaFold 结构中提取的残基一致。
- 移除具有内部冲突标签的突变。
- 对重复的突变条目进行去重。

### 最终数据集规模

文件：`data/processed/brca1_brct_clinvar_clean.csv`

| 划分 | 数量 |
| --- | ---: |
| 总变异数 | 201 |
| 致病 / 疑似致病 | 132 |
| 良性 / 疑似良性 | 69 |
| pLDDT < 70 的变异数 | 10 |
| pLDDT < 70 的变异比例 | 5.0% |
| pLDDT < 90 的变异数 | 38 |
| pLDDT < 90 的变异比例 | 18.9% |

后续阈值扫描中还生成了 `pLDDT < 90` 掩蔽版本，对应文件为 `data/processed/brca1_brct_plddt90_clinvar_clean.csv` 和 `data/processed/brca1_brct_plddt90_saprot_sequences.json`。两个版本使用同一批 201 个 BRCA1 BRCT ClinVar 错义突变，仅结构 token 掩蔽阈值发生变化。

## 3. 实验 1：零样本基线 (Zero-shot Baseline)

### 目的

评估在 BRCA1 BRCT 功能结构域内，与仅使用序列的蛋白质语言模型相比，加入结构 token 是否能提升零样本致病性排序表现。

### 对比模型

- `ESM-2`：仅基于序列的掩码语言模型（masked language model）。
- `SaProt full`：使用 BRCT 区域完整野生型 `Seq + 3Di` token 的结构感知模型。
- `SaProt masked`：根据 pLDDT，将 BRCT 区域中低于指定阈值的 `3Di` token 替换为 `#` 的结构感知模型。

### 评分定义

对于突变 `WT -> Mut`，零样本得分公式为：

```text
score = log P(WT | context) - log P(Mut | context)
```

因此，得分越大表明模型越强烈地倾向于野生型残基而非突变型残基，这被解释为该突变具有破坏性（deleterious）的更强证据。

### 主要结果

默认 pLDDT 掩蔽阈值：`70`。

文件：

- `results/brca1_brct_zero_shot_scores.csv`
- `results/brca1_brct_zero_shot_metrics.json`
- `results/brca1_brct_zero_shot_scores_roc.png`

| 模型 | ROC-AUC | ROC-AUC 95% CI | 平均精度 (AP) | AP 95% CI | 准确率 | F1 分数 |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| ESM-2 | 0.9168 | [0.8749, 0.9519] | 0.9563 | [0.9285, 0.9784] | 0.7114 | 0.8176 |
| SaProt full 3Di (完整) | 0.9223 | [0.8778, 0.9583] | 0.9562 | [0.9187, 0.9823] | 0.6965 | 0.8076 |
| SaProt masked 3Di, pLDDT < 70 | 0.9192 | [0.8722, 0.9582] | 0.9499 | [0.9013, 0.9825] | 0.6965 | 0.8063 |

### 简单偏置基线

为检验模型结果是否仅由结构置信度或区域分布偏置解释，额外加入两个不使用语言模型的简单 baseline：

- `pLDDT-only`：直接使用突变位点 pLDDT 作为致病性排序分数。
- `region-only`：在 BRCT domain-local 实验中没有有效区域差异，因此等价于常数基线。

| 基线 | ROC-AUC | ROC-AUC 95% CI | 平均精度 (AP) | AP 95% CI | 解释 |
| --- | ---: | ---: | ---: | ---: | --- |
| pLDDT-only | 0.7768 | [0.7037, 0.8426] | 0.8538 | [0.7878, 0.9133] | pLDDT 本身具有一定预测偏置，但弱于模型得分。 |
| region-only | 0.5000 | [0.5000, 0.5000] | 0.6567 | [0.5920, 0.7164] | 单一区域内的常数基线，不能解释模型性能。 |

### Bootstrap 差值检验

| 对比 | ROC-AUC 差值 | ROC-AUC 差值 95% CI | AP 差值 | AP 差值 95% CI |
| --- | ---: | ---: | ---: | ---: |
| SaProt full 3Di - ESM-2 | +0.0055 | [-0.0375, +0.0444] | -0.0001 | [-0.0379, +0.0307] |
| SaProt masked 3Di (70) - SaProt full 3Di | -0.0031 | [-0.0145, +0.0055] | -0.0063 | [-0.0197, +0.0021] |
| SaProt masked 3Di (70) - ESM-2 | +0.0024 | [-0.0412, +0.0419] | -0.0064 | [-0.0545, +0.0298] |

### 解释

BRCA1 BRCT 区域内，`SaProt full 3Di` 的 ROC-AUC 略高于 `ESM-2`，但 bootstrap 差值区间跨 0，因此不能声称结构 token 带来显著提升。默认阈值 `70` 下，masked 3Di 也没有优于完整 3Di。因此，BRCA1-BRCT 的阈值 `70` 结果应被描述为中性或轻微负向，而非支持性结果。

## 4. 实验 2：基于 pLDDT 的 3Di 掩蔽消融实验 (Ablation)

### 目的

测试在 BRCA1 BRCT 结构域局部上下文中，将中低置信度 AlphaFold 衍生结构 token 替换为未知结构 token，是否能改善结构条件突变评分。

### 阈值 70 下的主要结果

| 对比 | ROC-AUC | ROC-AUC 95% CI | 平均精度 (AP) | AP 95% CI |
| --- | ---: | ---: | ---: | ---: |
| SaProt full 3Di | 0.9223 | [0.8778, 0.9583] | 0.9562 | [0.9187, 0.9823] |
| SaProt masked 3Di, pLDDT < 70 | 0.9192 | [0.8722, 0.9582] | 0.9499 | [0.9013, 0.9825] |
| 掩蔽 - 完整 | -0.0031 | [-0.0145, +0.0055] | -0.0063 | [-0.0197, +0.0021] |

### 阈值 90 下的主要结果

文件：

- `results/brca1_brct_plddt90_zero_shot_scores.csv`
- `results/brca1_brct_plddt90_zero_shot_metrics.json`
- `results/brca1_brct_plddt90_zero_shot_scores_roc.png`

| 对比 | ROC-AUC | ROC-AUC 95% CI | 平均精度 (AP) | AP 95% CI |
| --- | ---: | ---: | ---: | ---: |
| ESM-2 | 0.9168 | [0.8749, 0.9519] | 0.9563 | [0.9285, 0.9784] |
| SaProt full 3Di | 0.9223 | [0.8778, 0.9583] | 0.9562 | [0.9187, 0.9823] |
| SaProt masked 3Di, pLDDT < 90 | 0.9300 | [0.8923, 0.9611] | 0.9693 | [0.9516, 0.9837] |

阈值 `90` 的 bootstrap 差值结果如下：

| 对比 | ROC-AUC 差值 | ROC-AUC 差值 95% CI | AP 差值 | AP 差值 95% CI |
| --- | ---: | ---: | ---: | ---: |
| SaProt full 3Di - ESM-2 | +0.0055 | [-0.0375, +0.0444] | -0.0001 | [-0.0379, +0.0307] |
| SaProt masked 3Di (90) - SaProt full 3Di | +0.0077 | [-0.0172, +0.0388] | +0.0131 | [-0.0067, +0.0460] |
| SaProt masked 3Di (90) - ESM-2 | +0.0132 | [-0.0162, +0.0439] | +0.0130 | [-0.0038, +0.0347] |

### 解释

阈值 `90` 下，BRCA1-BRCT 的 masked SaProt 取得最高 ROC-AUC 和 AP，方向上与 TP53/MSH2 的最佳阈值一致。然而，masked90 相对完整 3Di 和 ESM-2 的差值置信区间均跨 0，因此该结果只能被描述为正向趋势，而不是统计显著提升。BRCT 区域中 `pLDDT < 90` 的变异比例为 18.9%，低于 MSH2 的 46.0%，这可能限制了 masking 效果的可观察空间。

## 5. 区域级 Bootstrap 分析

### 目的

分析 BRCA1 BRCT 内部不同子区域是否存在不同的掩蔽收益。

### 区域划分

当前使用如下 BRCT 子区域划分：

- BRCT1：位置 1646-1736
- Linker：位置 1737-1755
- BRCT2：位置 1756-1859

### 阈值 70 的区域结果

文件：

- `results/brca1_brct_regional_bootstrap.csv`
- `results/brca1_brct_regional_bootstrap.json`

| 区域 | N | 致病 | 良性 | 平均 pLDDT | pLDDT < 70 | Full AUC | Masked AUC | Masked - Full AUC | AUC 差值 95% CI | Masked - Full AP | AP 差值 95% CI |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| BRCT1 1646-1736 | 93 | 59 | 34 | 91.9 | 7 | 0.9482 | 0.9531 | +0.0050 | [-0.0014, +0.0143] | +0.0027 | [-0.0003, +0.0074] |
| Linker 1737-1755 | 24 | 21 | 3 | 95.2 | 0 | 0.9683 | 0.9683 | +0.0000 | [0.0000, 0.0000] | +0.0000 | [-0.0000, +0.0000] |
| BRCT2 1756-1859 | 84 | 52 | 32 | 90.1 | 3 | 0.8816 | 0.8690 | -0.0126 | [-0.0420, +0.0072] | -0.0204 | [-0.0510, +0.0031] |

### 阈值 90 的区域结果

文件：

- `results/brca1_brct_plddt90_regional_bootstrap.csv`
- `results/brca1_brct_plddt90_regional_bootstrap.json`

| 区域 | N | 致病 | 良性 | 平均 pLDDT | pLDDT < 90 | Full AUC | Masked AUC | Masked - Full AUC | AUC 差值 95% CI | Masked - Full AP | AP 差值 95% CI |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| BRCT1 1646-1736 | 93 | 59 | 34 | 91.9 | 9 | 0.9482 | 0.9551 | +0.0070 | [-0.0112, +0.0253] | +0.0050 | [-0.0033, +0.0150] |
| Linker 1737-1755 | 24 | 21 | 3 | 95.2 | 1 | 0.9683 | 0.9048 | -0.0635 | [-0.2361, +0.0000] | -0.0103 | [-0.0427, +0.0000] |
| BRCT2 1756-1859 | 84 | 52 | 32 | 90.1 | 28 | 0.8816 | 0.9038 | +0.0222 | [-0.0300, +0.0834] | +0.0300 | [-0.0138, +0.0979] |

### 解释

BRCT 子区域分析没有发现统计显著的局部提升。阈值 `90` 下，BRCT1 和 BRCT2 都呈正向趋势，其中 BRCT2 的 AUC 从 `0.8816` 提升到 `0.9038`，AP 从 `0.9227` 提升到 `0.9526`，但差值区间跨 0。Linker 区域样本较少且良性样本只有 3 个，不适合作为稳定区域结论来源。

## 6. 实验 3：空间分割线性探测 (Spatial Split Linear Probe)

### 目的

在空间分割（而非随机分割）下提供轻量级的有监督验证，以减少来自相邻或空间接近残基位置的简单信息泄露。该实验不作为主要机制证据，而是检验 zero-shot 分数和 pLDDT 是否能在一个简单线性模型中提供可迁移的判别信号。

BRCA1-BRCT 线性探测使用阈值 `90` 的 zero-shot 结果，因为该阈值在 BRCT 实验中取得最高 AUC/AP。

### 特征

线性探测使用了以下特征：

- `esm2_score`
- `saprot_full_score`
- `saprot_masked_score`
- `plddt`

### 结果

文件：

- `results/brca1_brct_plddt90_linear_probe_metrics.json`
- `results/brca1_brct_plddt90_linear_probe_predictions.csv`

| 划分 | ROC-AUC | 平均精度 (AP) | 准确率 | F1 分数 |
| --- | ---: | ---: | ---: | ---: |
| 训练集 (Train) | 0.9155 | 0.9594 | 0.8500 | 0.8667 |
| 空间测试集 (Spatial test) | 0.9509 | 0.9782 | 0.8430 | 0.8914 |

### 系数

| 特征 | 系数 |
| --- | ---: |
| `esm2_score` | +0.4090 |
| `saprot_full_score` | -0.0583 |
| `saprot_masked_score` | +0.6116 |
| `plddt` | -0.0455 |

### 解释

BRCA1-BRCT 的空间分割线性探测在测试空间 cluster 上取得较高 ROC-AUC (`0.9509`) 和 AP (`0.9782`)，说明 zero-shot 分数在该 domain-local 设置下具有较强判别信号。最大的正系数出现在 `saprot_masked_score` 上，与 threshold `90` masked SaProt 在 zero-shot 中取得最高 AUC/AP 的趋势一致。

需要注意的是，线性探测是轻量级监督验证，不能替代 zero-shot bootstrap 结论。由于 BRCA1-BRCT 的 masked90 相对 full/ESM-2 的 zero-shot 差值区间跨 0，该线性探测结果应作为支持性信号，而不是主要统计证据。

## 7. 当前证据与项目主张的对照

### 已获支持

- BRCA1-BRCT domain-local scoring 已成功实现，避免了 BRCA1 全长超过模型输入长度的问题。
- BRCT 区域清洗后包含 201 个 ClinVar 错义突变，其中致病 132 个、良性 69 个。
- 阈值 `90` 下，`SaProt masked 3Di` 在 BRCT 区域取得最高 ROC-AUC 和 AP。
- `pLDDT-only` 和 `region-only` baseline 弱于模型得分，不能解释模型整体性能。
- 空间分割线性探测成功运行，并显示 `saprot_masked_score` 是线性模型中最大的正向特征。

### 部分获支持

- `SaProt full 3Di` 相对 `ESM-2` 只有很小正向 AUC 差值，且 bootstrap 差值区间跨 0。
- 阈值 `90` 的 masked SaProt 相对 full 和 ESM-2 呈正向趋势，但 bootstrap 差值区间跨 0，因此不能表述为显著提升。
- BRCT 子区域分析没有显著局部提升，仅 BRCT1/BRCT2 出现方向性正向趋势。

### 尚未完成

- BRCA1 其他功能区域，如 RING domain，尚未分析。
- BRCA1 的窗口化全长突变评分尚未实现。
- BRCA1 的 PyMOL 或 AlphaFold pLDDT 可视化尚未实现。

## 8. 建议的报告用语

当前 BRCA1-BRCT 结果应作为长蛋白功能结构域扩展实验，而不是强统计泛化证据：

```text
由于 BRCA1 全长超过当前 ESM-2/SaProt 输入长度限制，我们对其 C 端 BRCT repeats 区域 1646-1859 进行了 domain-local 分析。清洗后得到 201 个 ClinVar 错义突变，其中致病/疑似致病 132 个、良性/疑似良性 69 个。结果显示，pLDDT < 90 的 masked SaProt 在该区域取得最高 ROC-AUC 和 AP（AUC 0.9300，AP 0.9693），高于 SaProt full（AUC 0.9223，AP 0.9562）和 ESM-2（AUC 0.9168，AP 0.9563）。不过，pairwise bootstrap 差值区间仍跨 0，因此该结果应被描述为正向趋势而非显著提升。空间分割线性探测进一步显示，saprot_masked_score 是线性模型中最大的正向特征，并在测试空间 cluster 上取得较高 ROC-AUC 和 AP。总体而言，BRCA1-BRCT 实验说明该流程可以扩展到长蛋白的功能结构域局部分析，但当前最强统计支持仍来自 TP53 和 MSH2。
```
