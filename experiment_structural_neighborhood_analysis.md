# 补充实验四：机制分析 4c - 结构邻域分析

## 实验目的

本实验分析 pLDDT-aware masking 的收益是否与突变位点在三维结构中靠近低置信度区域有关。

原始假设是：

```text
如果一个突变位点在物理空间上更靠近低 pLDDT 区域，那么该突变位点更容易受到低置信度 3Di token 噪声干扰，因此 pLDDT-aware masking 的收益可能更大。
```

需要注意：单个突变没有 AUC，AUC 只能在一组样本上计算。因此本实验使用两层指标：

| 层面 | 指标 | 作用 |
|---|---|---|
| 单突变层面 | correct-direction score gain | 用于散点图，观察每个突变的 masking 收益与距离的关系 |
| 分箱层面 | binned AUC gain, 即 `masked AUC - full AUC` | 按距离分组后计算真正的 AUC 提升 |

## 方法

### 1. 最近低 pLDDT 区域距离

从 AlphaFold PDB 中提取每个残基 C-alpha 坐标和 pLDDT。对于每个突变位点，计算其到所有 `pLDDT < 90` 残基的最小三维欧氏距离：

```text
d_i = min_j || coord_i - coord_j ||_2, where pLDDT_j < 90
```

其中 `d_i = 0` 表示该突变位点自身就在低/中置信度区域中。

### 2. 单点 masking 收益

SaProt zero-shot score 越大，表示越倾向于致病。因此对于每个突变，定义 correct-direction score gain：

| 标签 | correct-direction score gain |
|---|---|
| pathogenic | `saprot_masked_score - saprot_full_score` |
| benign | `saprot_full_score - saprot_masked_score` |

该指标越大，表示 masking 后分数朝正确分类方向移动得越多。

### 3. 分箱 AUC 收益

按最近低 pLDDT 距离进行 quantile 分箱，并在每个箱内计算：

```text
AUC gain = AUC(saprot_masked_score) - AUC(saprot_full_score)
AP gain = AP(saprot_masked_score) - AP(saprot_full_score)
```

由于 TP53/MSH2 中有大量突变点自身位于低 pLDDT 区域，`qcut` 自动合并重复边界后实际得到 3 个距离箱。

## 实现与输出

新增脚本：

```text
scripts/run_structural_neighborhood_analysis.py
```

运行命令：

```bash
PYTHONPATH=. python scripts/run_structural_neighborhood_analysis.py --config configs/tp53_plddt90.yaml --output-prefix tp53_plddt90
PYTHONPATH=. python scripts/run_structural_neighborhood_analysis.py --config configs/msh2_plddt90.yaml --output-prefix msh2_plddt90
```

输出文件：

| 文件 | 内容 |
|---|---|
| `results/tp53_plddt90_structural_neighborhood.csv` | TP53 每个突变的最近低 pLDDT 距离和 score gain |
| `results/tp53_plddt90_structural_neighborhood_bins.csv` | TP53 距离分箱 AUC/AP gain |
| `results/tp53_plddt90_structural_neighborhood_metrics.json` | TP53 聚合相关性指标 |
| `results/tp53_plddt90_structural_neighborhood_scatter.png` | TP53 距离 vs 单点 score gain 散点图 |
| `results/tp53_plddt90_structural_neighborhood_binned_auc.png` | TP53 距离分箱 AUC gain 图 |
| `results/msh2_plddt90_structural_neighborhood.csv` | MSH2 每个突变的最近低 pLDDT 距离和 score gain |
| `results/msh2_plddt90_structural_neighborhood_bins.csv` | MSH2 距离分箱 AUC/AP gain |
| `results/msh2_plddt90_structural_neighborhood_metrics.json` | MSH2 聚合相关性指标 |
| `results/msh2_plddt90_structural_neighborhood_scatter.png` | MSH2 距离 vs 单点 score gain 散点图 |
| `results/msh2_plddt90_structural_neighborhood_binned_auc.png` | MSH2 距离分箱 AUC gain 图 |

## TP53 结果

TP53 数据集包含 373 个突变，其中 pathogenic/likely pathogenic 为 228 个，benign/likely benign 为 145 个。

### 聚合指标

| 指标 | 数值 |
|---|---:|
| mean nearest low-pLDDT distance | 9.7378 A |
| median nearest low-pLDDT distance | 10.5154 A |
| mean correct-direction score gain | +0.2170 |
| Spearman(distance, score gain) | +0.2296 |
| p-value | 7.52e-06 |
| Spearman(bin distance, AUC gain) | -0.5000 |
| p-value | 0.6667 |

### 距离分箱结果

| Distance bin | N | Pathogenic | Benign | Mean distance | Mean score gain | Full AUC | Masked AUC | AUC gain | AP gain |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| <= 10.515 A | 188 | 73 | 115 | 3.7531 | -0.0030 | 0.9521 | 0.9726 | +0.0205 | +0.0198 |
| 10.515-15.543 A | 94 | 80 | 14 | 12.6985 | +0.4838 | 0.9688 | 0.9607 | -0.0080 | -0.0014 |
| 15.543-26.915 A | 91 | 75 | 16 | 19.0437 | +0.3960 | 0.9333 | 0.9458 | +0.0125 | +0.0020 |

### TP53 解读

TP53 的结果不完全支持“距离越近，单点 score gain 越大”的简单假设。单点层面，最近低 pLDDT 距离与 correct-direction score gain 的 Spearman 相关为正，即距离越远，score gain 反而越大。

但从分箱 AUC 看，最近距离箱的 AUC gain 最大，为 `+0.0205`，说明靠近低 pLDDT 区域的样本组确实能从 masking 获得较明显的排序性能提升。由于只有 3 个有效距离箱，分箱相关性不显著，不能作强统计结论。

因此 TP53 应保守表述为：

```text
TP53 的分箱 AUC 结果与结构邻域假设部分一致，但单突变 score gain 相关性方向相反，说明 masking 收益还受到区域标签分布、突变类型和结构域背景等因素影响。
```

## MSH2 结果

MSH2 数据集包含 426 个突变，其中 pathogenic/likely pathogenic 为 141 个，benign/likely benign 为 285 个。

### 聚合指标

| 指标 | 数值 |
|---|---:|
| mean nearest low-pLDDT distance | 4.7494 A |
| median nearest low-pLDDT distance | 3.8377 A |
| mean correct-direction score gain | +0.5132 |
| Spearman(distance, score gain) | -0.2475 |
| p-value | 2.29e-07 |
| Spearman(bin distance, AUC gain) | -0.5000 |
| p-value | 0.6667 |

### 距离分箱结果

| Distance bin | N | Pathogenic | Benign | Mean distance | Mean score gain | Full AUC | Masked AUC | AUC gain | AP gain |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| <= 3.838 A | 215 | 40 | 175 | 0.3373 | +0.9564 | 0.7933 | 0.8434 | +0.0501 | +0.0687 |
| 3.838-8.796 A | 104 | 47 | 57 | 6.0982 | +0.1286 | 0.9634 | 0.9672 | +0.0037 | +0.0051 |
| 8.796-19.577 A | 107 | 54 | 53 | 12.3040 | -0.0037 | 0.9560 | 0.9623 | +0.0063 | +0.0034 |

### MSH2 解读

MSH2 结果较好地支持结构邻域假设。单突变层面，距离与 correct-direction score gain 呈显著负相关，Spearman rho 为 `-0.2475`，p-value 为 `2.29e-07`，说明越靠近低 pLDDT 区域，masking 后分数越倾向于朝正确方向移动。

分箱结果也与该趋势一致：最近距离箱的 AUC gain 为 `+0.0501`，AP gain 为 `+0.0687`，明显高于中远距离箱。该结果说明，MSH2 中靠近低置信度结构区域的突变确实更受益于 pLDDT-aware masking。

## 总体结论

结构邻域分析给出了一个更细粒度的机制视角：masking 的收益与突变位点和低置信度结构区域的三维关系有关，但这种关系具有蛋白依赖性。

| 观察 | TP53 | MSH2 |
|---|---|---|
| 单点距离 vs score gain | 正相关，和预期相反 | 显著负相关，符合预期 |
| 最近距离箱 AUC gain | 最大，部分支持假设 | 最大，强支持假设 |
| 分箱 AUC 趋势显著性 | 不显著 | 不显著，bin 数少 |
| 结论强度 | 部分支持，需要保守 | 支持较强 |

因此，4c 的结论应写得保守：

```text
MSH2 上，突变位点越靠近低 pLDDT 区域，masking 收益越大，支持“低置信度结构邻域噪声”这一机制解释。TP53 上，最近距离箱的 AUC gain 也最高，但单点 score gain 与距离呈相反趋势，说明结构邻域不是唯一决定因素，masking 收益还受到蛋白区域、标签分布和突变背景影响。
```

## 与 4a/4b 的关系

| 机制实验 | 主要发现 |
|---|---|
| 4a Attention 分析 | masking 后注意力从低 pLDDT 区域转向高 pLDDT 区域 |
| 4b Embedding 分析 | masking 后 pathogenic/benign 的 embedding 几何分离改善 |
| 4c 结构邻域分析 | MSH2 中靠近低 pLDDT 区域的突变更受益；TP53 部分支持但不完全一致 |

三者合起来支持一个保守机制解释：

```text
pLDDT-aware masking 通过抑制低置信度结构 token 的影响，改变模型对结构上下文的注意力和突变位点表征；在某些蛋白中，尤其是 MSH2，靠近低置信度结构区域的突变从 masking 中获益更明显。
```

## 推荐报告表述

```text
We further examined whether the benefit of pLDDT-aware masking depends on the 3D structural proximity between a mutation site and low-confidence AlphaFold regions. For each variant, we computed the minimum C-alpha distance to residues with pLDDT < 90 and compared this distance with the masking-induced score change. On MSH2, distance was significantly negatively correlated with the correct-direction masking gain, and the nearest-distance bin showed the largest AUC/AP improvement. This supports the hypothesis that low-confidence structural neighborhoods can introduce noise that is mitigated by masking. On TP53, the nearest-distance bin also showed the largest AUC gain, but the single-variant score-gain correlation had the opposite direction, indicating that structural proximity is not the only factor controlling masking benefit. We therefore treat structural neighborhood effects as supportive but protein-dependent mechanistic evidence.
```
