# 实验总结：全长 MSH2 突变效应预测

本文档总结了将 `project.md` 中的零样本突变效应预测流程迁移到全长 MSH2 后得到的实验结果。当前实验范围是全长人类 MSH2，残基位置为 1-934。

## 1. 项目目标

本实验的目标是检验在 TP53 之外，结构感知（structure-aware）的蛋白质语言模型是否仍能提升零样本（zero-shot）突变效应预测能力，以及基于 pLDDT 掩蔽 AlphaFold 结构 token 的策略是否具有跨蛋白支持。

当前核心问题是：

- 在全长 MSH2 的零样本致病性排序中，`SaProt (Seq + 3Di)` 是否优于仅使用序列的 `ESM-2`？
- 基于 pLDDT 对中低置信度 `3Di` token 进行掩蔽，是否能使 `SaProt` 的表现优于完整未掩蔽的 `3Di`？
- 掩蔽收益是否集中在 MSH2 的特定区域，从而支持“低置信度结构 token 可能引入噪声”的机制解释？

## 2. 数据集与预处理

### 目标对象

- 蛋白质：全长人类 MSH2
- UniProt ID：`P43246`
- 残基范围：`1-934`
- 结构来源：AlphaFold 模型 `AF-P43246-F1-model_v6.pdb`
- 结构分词（Tokenization）：`Foldseek structureto3didescriptor`

### ClinVar 数据清洗规则

数据集基于 ClinVar 的 `variant_summary.txt.gz` 构建，并使用与 TP53 实验一致的过滤条件：

- 仅保留 `GeneSymbol == MSH2` 的记录。
- 仅解析单氨基酸替换（single amino-acid substitutions）。
- 移除 uncertain、conflicting、not-provided、association-only、drug-response、risk-factor、protective 以及其他非二分类的临床标签。
- 将 `Pathogenic/Likely pathogenic`（致病/疑似致病）合并为标签 `1`。
- 将 `Benign/Likely benign`（良性/疑似良性）合并为标签 `0`。
- 要求 ClinVar 中的野生型残基与 UniProt 参考序列一致。
- 要求野生型残基与从 AlphaFold 结构中提取的残基一致。
- 移除具有内部冲突标签的突变。
- 对重复的突变条目进行去重。

### 最终数据集规模

文件：`data/processed/msh2_clinvar_clean.csv`

| 划分 | 数量 |
| --- | ---: |
| 总变异数 | 426 |
| 致病 / 疑似致病 | 141 |
| 良性 / 疑似良性 | 285 |
| pLDDT < 70 的变异数 | 40 |
| pLDDT < 70 的变异比例 | 9.4% |
| pLDDT < 90 的变异数 | 196 |
| pLDDT < 90 的变异比例 | 46.0% |

MSH2 的标签分布明显优于 PTEN，因此更适合作为 TP53 之外的跨蛋白补充实验。后续还生成了 `pLDDT < 90` 掩蔽版本，对应文件为 `data/processed/msh2_plddt90_clinvar_clean.csv` 和 `data/processed/msh2_plddt90_saprot_sequences.json`。两个阈值版本使用同一批 426 个 MSH2 ClinVar 错义突变，仅结构 token 掩蔽阈值发生变化。

## 3. 实验 1：零样本基线 (Zero-shot Baseline)

### 目的

评估与仅使用序列的蛋白质语言模型相比，加入结构 token 是否能提升全长 MSH2 上的零样本致病性排序表现。

### 对比模型

- `ESM-2`：仅基于序列的掩码语言模型（masked language model）。
- `SaProt full`：使用完整野生型 `Seq + 3Di` token 的结构感知模型。
- `SaProt masked`：根据 pLDDT，将低于指定阈值的 `3Di` token 替换为 `#` 的结构感知模型。

### 评分定义

对于突变 `WT -> Mut`，零样本得分公式为：

```text
score = log P(WT | context) - log P(Mut | context)
```

因此，得分越大表明模型越强烈地倾向于野生型残基而非突变型残基，这被解释为该突变具有破坏性（deleterious）的更强证据。

### 主要结果

默认 pLDDT 掩蔽阈值：`70`。

文件：

- `results/msh2_zero_shot_scores.csv`
- `results/msh2_zero_shot_metrics.json`
- `results/msh2_zero_shot_scores_roc.png`

| 模型 | ROC-AUC | ROC-AUC 95% CI | 平均精度 (AP) | AP 95% CI | 准确率 | F1 分数 |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| ESM-2 | 0.8870 | [0.8482, 0.9216] | 0.7979 | [0.7219, 0.8665] | 0.3779 | 0.5138 |
| SaProt full 3Di (完整) | 0.9105 | [0.8747, 0.9421] | 0.8772 | [0.8193, 0.9228] | 0.4272 | 0.5361 |
| SaProt masked 3Di, pLDDT < 70 | 0.9070 | [0.8699, 0.9394] | 0.8752 | [0.8173, 0.9205] | 0.4225 | 0.5323 |

### 简单偏置基线

为检验模型结果是否仅由结构置信度或区域分布偏置解释，额外加入两个不使用语言模型的简单 baseline：

- `pLDDT-only`：直接使用突变位点 pLDDT 作为致病性排序分数。
- `region-only`：仅根据粗粒度位置区域打分。

| 基线 | ROC-AUC | ROC-AUC 95% CI | 平均精度 (AP) | AP 95% CI | 解释 |
| --- | ---: | ---: | ---: | ---: | --- |
| pLDDT-only | 0.7692 | [0.7228, 0.8147] | 0.6062 | [0.5185, 0.6957] | pLDDT 本身具有一定预测偏置，但明显弱于 SaProt。 |
| region-only | 0.4624 | [0.4204, 0.5022] | 0.3178 | [0.2750, 0.3641] | 粗粒度区域位置不能解释 SaProt 的性能。 |

### Bootstrap 差值检验

| 对比 | ROC-AUC 差值 | ROC-AUC 差值 95% CI | AP 差值 | AP 差值 95% CI |
| --- | ---: | ---: | ---: | ---: |
| SaProt full 3Di - ESM-2 | +0.0235 | [+0.0050, +0.0414] | +0.0793 | [+0.0350, +0.1276] |
| SaProt masked 3Di (70) - SaProt full 3Di | -0.0036 | [-0.0083, +0.0004] | -0.0020 | [-0.0101, +0.0056] |
| SaProt masked 3Di (70) - ESM-2 | +0.0200 | [+0.0009, +0.0383] | +0.0773 | [+0.0337, +0.1252] |

### 解释

全长 MSH2 零样本基线支持结构感知模型的跨蛋白有效性：`SaProt full 3Di` 在 ROC-AUC 和 AP 上显著优于 `ESM-2`，且 bootstrap 差值区间均不跨 0。默认阈值 `70` 下，`SaProt masked` 仍显著优于 `ESM-2`，但相对完整 3Di 略低，且差值区间跨 0。因此，MSH2 上不应将阈值 `70` 描述为改善完整 3Di 的有效掩蔽阈值。

## 4. 实验 2：基于 pLDDT 的 3Di 掩蔽消融实验 (Ablation)

### 目的

测试将中低置信度的 AlphaFold 衍生结构 token 替换为未知结构 token，是否能减少结构感知突变评分中的噪声。

### 阈值 70 下的主要全长结果

| 对比 | ROC-AUC | ROC-AUC 95% CI | 平均精度 (AP) | AP 95% CI |
| --- | ---: | ---: | ---: | ---: |
| SaProt full 3Di | 0.9105 | [0.8747, 0.9421] | 0.8772 | [0.8193, 0.9228] |
| SaProt masked 3Di, pLDDT < 70 | 0.9070 | [0.8699, 0.9394] | 0.8752 | [0.8173, 0.9205] |
| 掩蔽 - 完整 | -0.0036 | [-0.0083, +0.0004] | -0.0020 | [-0.0101, +0.0056] |

### 阈值 90 下的主要全长结果

文件：

- `results/msh2_plddt90_zero_shot_scores.csv`
- `results/msh2_plddt90_zero_shot_metrics.json`
- `results/msh2_plddt90_zero_shot_scores_roc.png`

| 对比 | ROC-AUC | ROC-AUC 95% CI | 平均精度 (AP) | AP 95% CI |
| --- | ---: | ---: | ---: | ---: |
| ESM-2 | 0.8870 | [0.8482, 0.9216] | 0.7979 | [0.7219, 0.8665] |
| SaProt full 3Di | 0.9105 | [0.8747, 0.9421] | 0.8772 | [0.8193, 0.9228] |
| SaProt masked 3Di, pLDDT < 90 | 0.9278 | [0.8955, 0.9549] | 0.9125 | [0.8734, 0.9441] |

阈值 `90` 的 bootstrap 差值结果如下：

| 对比 | ROC-AUC 差值 | ROC-AUC 差值 95% CI | AP 差值 | AP 差值 95% CI |
| --- | ---: | ---: | ---: | ---: |
| SaProt full 3Di - ESM-2 | +0.0235 | [+0.0050, +0.0414] | +0.0793 | [+0.0350, +0.1276] |
| SaProt masked 3Di (90) - SaProt full 3Di | +0.0173 | [+0.0047, +0.0298] | +0.0354 | [+0.0068, +0.0715] |
| SaProt masked 3Di (90) - ESM-2 | +0.0408 | [+0.0209, +0.0611] | +0.1146 | [+0.0611, +0.1704] |

### 解释

MSH2 上阈值 `90` 是强阳性结果。`SaProt masked 3Di (90)` 同时显著优于 `SaProt full 3Di` 和 `ESM-2`，并且 ROC-AUC 与 AP 的 bootstrap 差值区间均不跨 0。这与 TP53 中阈值 `90` 最强的观察一致，支持如下假设：不仅低置信度结构 token，部分中等置信度 AlphaFold 结构 token 也可能在结构条件突变评分中引入噪声，而 pLDDT-aware masking 可以缓解这一问题。

## 5. 区域级 Bootstrap 分析

### 目的

分析 MSH2 上掩蔽收益来自哪些区域，并检验全局提升是否由特定结构区域驱动。

### 区域划分

当前使用三段粗粒度位置划分：

- N 端区域：位置 1-300
- 中部区域：位置 301-650
- C 端区域：位置 651-934

该划分用于统计分析，不等同于严格的结构域注释。后续若用于正式论文或报告，可再替换为基于文献的 MSH2 结构域边界。

### 阈值 70 的区域结果

文件：

- `results/msh2_regional_bootstrap.csv`
- `results/msh2_regional_bootstrap.json`

| 区域 | N | 致病 | 良性 | 平均 pLDDT | pLDDT < 70 | Full AUC | Masked AUC | Masked - Full AUC | AUC 差值 95% CI | Masked - Full AP | AP 差值 95% CI |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| N 端 1-300 | 143 | 28 | 115 | 82.9 | 22 | 0.9245 | 0.9245 | +0.0000 | [-0.0130, +0.0119] | +0.0021 | [-0.0480, +0.0532] |
| 中部 301-650 | 172 | 64 | 108 | 91.7 | 5 | 0.8963 | 0.8919 | -0.0043 | [-0.0107, +0.0003] | -0.0016 | [-0.0061, +0.0031] |
| C 端 651-934 | 111 | 49 | 62 | 85.3 | 13 | 0.9134 | 0.9115 | -0.0020 | [-0.0120, +0.0078] | -0.0010 | [-0.0176, +0.0148] |

### 阈值 90 的区域结果

文件：

- `results/msh2_plddt90_regional_bootstrap.csv`
- `results/msh2_plddt90_regional_bootstrap.json`

| 区域 | N | 致病 | 良性 | 平均 pLDDT | pLDDT < 90 | Full AUC | Masked AUC | Masked - Full AUC | AUC 差值 95% CI | Masked - Full AP | AP 差值 95% CI |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| N 端 1-300 | 143 | 28 | 115 | 82.9 | 92 | 0.9245 | 0.9335 | +0.0090 | [-0.0151, +0.0344] | +0.0539 | [-0.0084, +0.1486] |
| 中部 301-650 | 172 | 64 | 108 | 91.7 | 47 | 0.8963 | 0.9036 | +0.0074 | [-0.0094, +0.0248] | +0.0044 | [-0.0145, +0.0230] |
| C 端 651-934 | 111 | 49 | 62 | 85.3 | 57 | 0.9134 | 0.9437 | +0.0303 | [+0.0045, +0.0622] | +0.0570 | [+0.0059, +0.1318] |

### 解释

MSH2 的区域级 bootstrap 显示，阈值 `90` 的全局提升主要由 C 端区域 `651-934` 驱动。该区域中，masked 3Di 的 ROC-AUC 从 `0.9134` 提升到 `0.9437`，AUC 差值为 `+0.0303`，95% CI 为 `[+0.0045, +0.0622]`；AP 从 `0.8984` 提升到 `0.9554`，AP 差值为 `+0.0570`，95% CI 为 `[+0.0059, +0.1318]`。N 端和中部区域在阈值 `90` 下也呈正向趋势，但差值区间跨 0，因此不能表述为显著改善。

阈值 `70` 下，各区域均未显示稳定改善。这说明 MSH2 上的有效掩蔽并不是简单地去除最低置信度结构 token，而更可能需要掩蔽一部分中等置信度 token。

## 6. 实验 3：空间分割线性探测 (Spatial Split Linear Probe)

### 目的

在空间分割（而非随机分割）下提供轻量级的有监督验证，以减少来自相邻或空间接近残基位置的简单信息泄露。该实验不作为主要机制证据，而是检验 zero-shot 分数和 pLDDT 是否能在一个简单线性模型中提供可迁移的判别信号。

MSH2 线性探测使用阈值 `90` 的 zero-shot 结果，因为该阈值在 MSH2 的主实验中表现最好。

### 特征

线性探测使用了以下特征：

- `esm2_score`
- `saprot_full_score`
- `saprot_masked_score`
- `plddt`

### 结果

文件：

- `results/msh2_plddt90_linear_probe_metrics.json`
- `results/msh2_plddt90_linear_probe_predictions.csv`

| 划分 | ROC-AUC | 平均精度 (AP) | 准确率 | F1 分数 |
| --- | ---: | ---: | ---: | ---: |
| 训练集 (Train) | 0.9061 | 0.8543 | 0.8363 | 0.7087 |
| 空间测试集 (Spatial test) | 0.9486 | 0.9495 | 0.8450 | 0.8306 |

### 系数

| 特征 | 系数 |
| --- | ---: |
| `esm2_score` | +0.0748 |
| `saprot_full_score` | -0.1976 |
| `saprot_masked_score` | +0.6292 |
| `plddt` | -0.0002 |

### 解释

MSH2 的空间分割线性探测在测试空间 cluster 上取得了较高的 ROC-AUC (`0.9486`) 和 AP (`0.9495`)，说明 zero-shot 分数在空间分割设置下仍包含可迁移的判别信号。最大的正系数出现在 `saprot_masked_score` 上，这与阈值 `90` 的 masked SaProt 在零样本实验中表现最强相一致。

需要注意的是，训练集 ROC-AUC 低于测试集 ROC-AUC，说明该空间划分下两个 cluster 的难度并不完全一致。因此，该实验应被描述为轻量级监督验证，而不是完整的端到端监督微调结论。

## 7. 当前证据与项目主张的对照

### 已获支持

- 全长 MSH2 零样本实验已完成。
- MSH2 数据集标签分布较好，包含 426 个变异，其中致病 141 个、良性 285 个。
- 在 ROC-AUC 和 AP 方面，`SaProt full 3Di` 显著优于 `ESM-2`，且 bootstrap 差值区间不跨 0。
- 简单的 `pLDDT-only` 和 `region-only` baseline 不能解释 SaProt 的性能优势。
- 阈值 `90` 下，`SaProt masked 3Di` 显著优于 `SaProt full 3Di` 和 `ESM-2`。
- 区域级 bootstrap 显示，阈值 `90` 的改进主要由 MSH2 C 端区域 `651-934` 驱动，且该区域 AUC/AP 差值区间均不跨 0。
- 空间分割线性探测成功运行，并显示 `saprot_masked_score` 是线性模型中最大的正向特征。

### 部分获支持

- 阈值 `70` 下，`SaProt masked 3Di` 仍显著优于 `ESM-2`，但相对 `SaProt full 3Di` 没有改善，因此默认阈值 `70` 不应被表述为 MSH2 上的有效提升阈值。
- 当前区域划分是粗粒度位置划分，不是严格结构域划分；若用于正式展示，建议进一步依据 MSH2 文献结构域边界修订。

### 尚未完成

- MSH2 的更细粒度结构域边界分析尚未实现。
- MSH2 的 PyMOL 或 AlphaFold pLDDT 可视化尚未实现。

## 8. 建议的报告用语

当前 MSH2 结果可以作为 TP53 之外的跨蛋白补充证据，但仍应保守报告：

```text
为了检验 TP53 上观察到的现象是否可以迁移到其他蛋白，我们将同一零样本流程应用于全长 MSH2。清洗后的 MSH2 ClinVar 数据包含 426 个错义突变，其中致病/疑似致病 141 个，良性/疑似良性 285 个，标签分布明显优于 PTEN。结果显示，SaProt full 3Di 显著优于仅使用序列的 ESM-2，说明结构条件突变评分在 TP53 之外仍然有效。更重要的是，在 pLDDT < 90 掩蔽阈值下，masked SaProt 进一步显著优于完整 3Di 输入，ROC-AUC 从 0.9105 提升到 0.9278，AP 从 0.8772 提升到 0.9125，且 bootstrap 差值区间均不跨 0。区域级 bootstrap 显示，该提升主要来自 MSH2 C 端区域 651-934，其中 ROC-AUC 从 0.9134 提升到 0.9437，AP 从 0.8984 提升到 0.9554。相比之下，阈值 70 并未带来稳定改善。这一跨蛋白结果支持我们的核心假设：AlphaFold 结构 token 中的中低置信度信息可能成为结构条件突变评分的噪声来源，而基于 pLDDT 的结构 token 掩蔽可以缓解这一问题。
```

不应声称 MSH2 已经证明该策略在所有蛋白上普遍有效。更稳妥的定位是：MSH2 提供了一个标签分布较好的跨蛋白支持案例，并且独立复现了 `pLDDT < 90` 掩蔽优于完整 3Di 的现象。
