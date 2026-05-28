# 补充实验四：机制分析 4a - Attention 分析

## 实验目的

前面的实验已经表明，pLDDT-aware masking 能够提升 TP53 和 MSH2 上的 SaProt zero-shot 突变效应预测性能，并且这种提升不是由随机 token 删除造成的。本实验进一步分析模型内部注意力是否发生了符合机制假设的变化。

核心问题是：

```text
在 pLDDT-aware masking 后，SaProt 是否减少了对低置信度结构区域的注意力，并相对更多关注高置信度结构区域或突变位点附近区域？
```

需要强调的是，attention 分析只能作为机制线索，不能单独证明因果机制。因此本实验采用定量聚合统计加少量 case heatmap 的方式，避免只依赖主观挑选某个 head 或某个样本。

## 实验设计

分析对象仍然选择两个主力蛋白：

| 蛋白 | 作用 |
|---|---|
| TP53 | 主实验蛋白，pLDDT-aware masking 在 threshold 90 下显著提升 |
| MSH2 | 跨蛋白全长验证，pLDDT-aware masking 在 threshold 90 下显著提升 |

对于每个蛋白，从 zero-shot 结果中抽取 80 个代表性变异：

```text
按 masked - full 的 correct-direction score change 排序，并保持 pathogenic/benign 类别平衡。
```

其中 correct-direction change 的定义为：

| 标签 | correct-direction change |
|---|---|
| pathogenic | `saprot_masked_score - saprot_full_score` 越大越好 |
| benign | `saprot_full_score - saprot_masked_score` 越大越好 |

对于每个变异，使用与 zero-shot scoring 一致的输入：

```text
mask 掉突变位点 amino-acid token，然后提取该突变位点作为 query 时对所有 residue token 的 attention 分布。
```

attention 统计采用所有 layer 和所有 head 的平均值，而不是手动挑选特定 head。

## 评价指标

| 指标 | 含义 | 期望方向 |
|---|---|---|
| low pLDDT attention mass | 突变位点 query 分配到 pLDDT < 90 残基的 attention 总量 | masking 后下降 |
| high pLDDT attention mass | 突变位点 query 分配到 pLDDT >= 90 残基的 attention 总量 | masking 后上升 |
| local attention mass | 突变位点 +/- 10 个残基窗口内的 attention 总量 | 可能上升 |
| attention entropy | attention 分布熵，反映注意力是否更分散 | masking 后可能下降 |

其中 low/high pLDDT attention mass 是本实验最核心的机制指标。

## 实现与输出

新增脚本：

```text
scripts/run_attention_analysis.py
```

运行命令：

```bash
PYTHONPATH=. python scripts/run_attention_analysis.py --config configs/tp53_plddt90.yaml --output-prefix tp53_plddt90 --max-variants 80 --n-cases 2
PYTHONPATH=. python scripts/run_attention_analysis.py --config configs/msh2_plddt90.yaml --output-prefix msh2_plddt90 --max-variants 80 --n-cases 2
```

输出文件：

| 文件 | 内容 |
|---|---|
| `results/tp53_plddt90_attention_analysis.csv` | TP53 每个变异的 attention 统计 |
| `results/tp53_plddt90_attention_analysis_metrics.json` | TP53 聚合 attention 指标 |
| `results/tp53_plddt90_attention_case*.png` | TP53 代表性 case heatmap |
| `results/msh2_plddt90_attention_analysis.csv` | MSH2 每个变异的 attention 统计 |
| `results/msh2_plddt90_attention_analysis_metrics.json` | MSH2 聚合 attention 指标 |
| `results/msh2_plddt90_attention_case*.png` | MSH2 代表性 case heatmap |

## TP53 结果

TP53 抽取 80 个代表性变异进行 attention 分析。

| 指标 | SaProt full | SaProt pLDDT-mask | masked - full | 变化方向 |
|---|---:|---:|---:|---|
| low pLDDT attention mass | 0.3524 | 0.3149 | -0.0375 | 下降 |
| high pLDDT attention mass | 0.5974 | 0.6364 | +0.0389 | 上升 |
| local attention mass | 0.3426 | 0.3594 | +0.0168 | 上升 |
| attention entropy | 5.1885 | 5.1456 | -0.0430 | 下降 |

### TP53 解读

TP53 上，masking 后突变位点对低 pLDDT 残基的平均 attention mass 从 `0.3524` 降至 `0.3149`，同时对高 pLDDT 残基的 attention mass 从 `0.5974` 升至 `0.6364`。这与“masking 抑制低置信度结构噪声，使模型更多依赖高置信度结构信息”的机制假设一致。

此外，local attention mass 略有上升，attention entropy 下降，说明 masking 后突变位点的注意力分布略微更集中，并且更偏向局部结构上下文。

## MSH2 结果

MSH2 同样抽取 80 个代表性变异进行 attention 分析。

| 指标 | SaProt full | SaProt pLDDT-mask | masked - full | 变化方向 |
|---|---:|---:|---:|---|
| low pLDDT attention mass | 0.4447 | 0.4164 | -0.0283 | 下降 |
| high pLDDT attention mass | 0.5340 | 0.5616 | +0.0277 | 上升 |
| local attention mass | 0.3147 | 0.3228 | +0.0081 | 上升 |
| attention entropy | 5.5919 | 5.5776 | -0.0143 | 下降 |

### MSH2 解读

MSH2 上观察到与 TP53 一致的方向：pLDDT-aware masking 后，突变位点对低 pLDDT 区域的 attention mass 下降，对高 pLDDT 区域的 attention mass 上升。local attention mass 小幅上升，attention entropy 小幅下降。

由于 MSH2 是独立的全长蛋白验证集，该结果说明 attention 重分配现象不是 TP53 特例，而是在另一个蛋白上也能复现。

## 与 4b Embedding 分析的关系

Attention 分析和 embedding 分析从两个层面提供机制证据：

| 分析 | 观察 | 解释 |
|---|---|---|
| 4a Attention | masking 后低 pLDDT attention mass 下降，高 pLDDT attention mass 上升 | 模型内部注意力从不可靠结构区域转向更可靠结构区域 |
| 4b Embedding | masking 后 Fisher ratio、silhouette 等几何分离指标提高 | 突变位点表征更利于区分 pathogenic/benign |

二者合起来支持一个更完整的机制解释：

```text
pLDDT-aware masking 先改变模型对结构上下文的注意力分配，减少对低置信度结构 token 的依赖；这种注意力重分配进一步反映到突变位点 embedding 空间中，使 pathogenic 与 benign 的表征分离更清晰。
```

## 局限性

本实验需要保守解释，原因如下：

| 局限 | 说明 |
|---|---|
| attention 不等于因果解释 | attention weight 只能提供模型内部行为线索，不能单独证明因果机制 |
| 使用 layer/head 平均 | 平均化提高稳定性，但可能掩盖个别 head 的特殊功能 |
| 抽样分析而非全量分析 | 为控制计算量，每个蛋白抽取 80 个代表性变异 |
| case heatmap 仅作可视化 | 单个 heatmap 不能替代聚合统计 |

因此，attention 分析应作为辅助机制证据，而不是主性能证据。主性能结论仍应基于 zero-shot bootstrap、严格 masking 消融和 MSH2 跨蛋白验证。

## 总体结论

TP53 和 MSH2 的 attention 分析均显示，pLDDT-aware masking 后：

| 现象 | TP53 | MSH2 | 是否一致 |
|---|---:|---:|---|
| low pLDDT attention mass 下降 | -0.0375 | -0.0283 | 是 |
| high pLDDT attention mass 上升 | +0.0389 | +0.0277 | 是 |
| local attention mass 上升 | +0.0168 | +0.0081 | 是 |
| attention entropy 下降 | -0.0430 | -0.0143 | 是 |

这说明 masking 后，SaProt 的突变位点注意力确实从低置信度结构区域转向更高置信度和更局部的结构上下文。该结果为 pLDDT-aware masking 的有效性提供了机制层面的支持。

## 推荐报告表述

```text
To probe how pLDDT-aware masking changes model behavior, we analyzed SaProt attention maps before and after masking. For each selected variant, we used the mutation-site token as the query and averaged attention weights across all layers and heads. On both TP53 and MSH2, masking reduced the attention mass assigned to low-pLDDT residues and increased the attention mass assigned to high-pLDDT residues. Local attention around the mutation site also slightly increased, while attention entropy decreased. These observations suggest that pLDDT-aware masking redirects the model away from unreliable structural regions toward more reliable and locally relevant structural context. We treat this analysis as supportive mechanistic evidence rather than a standalone causal proof.
```
