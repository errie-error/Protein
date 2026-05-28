# 补充实验四：机制分析 4b - Embedding 空间分析

## 实验目的

前面的主实验和 masking 消融已经证明：pLDDT-aware masking 可以提升 TP53 和 MSH2 的零样本突变效应预测性能，并且这种提升不是由随机删除结构 token 或 token 数量减少造成的。

本实验进一步分析 masking 为什么有效。核心问题是：

```text
pLDDT-aware masking 是否改变了 SaProt 的突变表征空间，使 pathogenic 与 benign 变异更容易区分？
```

如果 masking 的作用机制是减少不可靠结构 token 带来的噪声，那么在 masked SaProt embedding 中，pathogenic 与 benign 两类突变应该表现出更好的几何分离，例如类间距离增大、类内方差减小、Fisher ratio 上升，或者轻量线性分类器更容易区分两类。

## 实验设计

本实验比较两种 SaProt 输入：

| 表征 | 描述 |
|---|---|
| SaProt full embedding | 使用完整 amino-acid + 3Di 序列，不进行 pLDDT-aware masking |
| SaProt pLDDT-mask embedding | 使用 pLDDT < 90 的结构 token masking 序列 |

为了与 zero-shot scoring 保持一致，每个突变的 embedding 定义为：

```text
将突变位点的 amino-acid token mask 后，提取 SaProt 最后一层 hidden state 中该突变位点对应 token 的向量。
```

这样 full 与 masked 的差异只来自结构 token masking，而不是来自不同的突变编码方式。

## 评价指标

| 指标 | 含义 | 期望方向 |
|---|---|---|
| centroid distance | pathogenic 与 benign 两类中心的欧氏距离 | 越大越好 |
| mean within-class variance | 两类样本到各自类别中心的平均平方距离 | 越小越好 |
| normalized centroid distance | 类间距离除以类内标准差 | 越大越好 |
| Fisher ratio | 类间距离平方 / 类内方差 | 越大越好 |
| silhouette score | 使用真实标签衡量 embedding 聚类分离度 | 越大越好 |
| linear probe CV ROC-AUC/AP | 在 embedding 上训练轻量 logistic regression 的 5-fold CV 表现 | 越大越好 |

需要注意：embedding 空间分析是机制证据，不是主性能指标。主性能结论仍以 zero-shot AUC/AP 和 bootstrap 差值为准。

## 实现与输出

新增脚本：

```text
scripts/run_embedding_analysis.py
```

运行命令：

```bash
PYTHONPATH=. python scripts/run_embedding_analysis.py --config configs/tp53_plddt90.yaml --output-prefix tp53_plddt90
PYTHONPATH=. python scripts/run_embedding_analysis.py --config configs/msh2_plddt90.yaml --output-prefix msh2_plddt90
```

输出文件：

| 文件 | 内容 |
|---|---|
| `results/tp53_plddt90_embedding_analysis_metrics.json` | TP53 embedding 指标 |
| `results/tp53_plddt90_embedding_tsne.png` | TP53 t-SNE 可视化 |
| `results/tp53_plddt90_embedding_analysis_summary.csv` | TP53 embedding summary 表 |
| `results/msh2_plddt90_embedding_analysis_metrics.json` | MSH2 embedding 指标 |
| `results/msh2_plddt90_embedding_tsne.png` | MSH2 t-SNE 可视化 |
| `results/msh2_plddt90_embedding_analysis_summary.csv` | MSH2 embedding summary 表 |

## TP53 结果

TP53 数据集包含 373 个 ClinVar missense variants，其中 pathogenic/likely pathogenic 为 228 个，benign/likely benign 为 145 个。

| 指标 | SaProt full | SaProt pLDDT-mask | masked - full | 变化方向 |
|---|---:|---:|---:|---|
| centroid distance | 2.2367 | 2.5425 | +0.3058 | 提升 |
| mean within-class variance | 16.8733 | 17.4947 | +0.6214 | 略升 |
| normalized centroid distance | 0.5445 | 0.6079 | +0.0634 | 提升 |
| Fisher ratio | 0.2965 | 0.3695 | +0.0730 | 提升 |
| silhouette score | 0.0533 | 0.0685 | +0.0151 | 提升 |
| linear probe CV ROC-AUC | 0.9457 | 0.9434 | -0.0023 | 基本持平，轻微下降 |
| linear probe CV AP | 0.9612 | 0.9589 | -0.0023 | 基本持平，轻微下降 |

### TP53 解读

TP53 上，pLDDT-mask 后 pathogenic 与 benign 的类别中心距离增大，normalized centroid distance、Fisher ratio 和 silhouette score 均上升。这说明从几何结构看，masked embedding 中两类突变的分离程度更强。

不过，TP53 的类内方差略有上升，linear probe CV AUC/AP 也轻微下降。因此 TP53 的 embedding 机制证据应当保守表述为：

```text
masking 改善了 TP53 embedding 空间中的类别几何分离，但没有带来线性 probe 指标的同步提升。
```

这与 TP53 zero-shot 主结果并不矛盾，因为 zero-shot scoring 使用的是 masked-LM 似然差，而不是直接在 embedding 上训练分类器。

## MSH2 结果

MSH2 数据集包含 426 个 ClinVar missense variants，其中 pathogenic/likely pathogenic 为 141 个，benign/likely benign 为 285 个。

| 指标 | SaProt full | SaProt pLDDT-mask | masked - full | 变化方向 |
|---|---:|---:|---:|---|
| centroid distance | 1.7315 | 1.9033 | +0.1718 | 提升 |
| mean within-class variance | 21.4992 | 18.9480 | -2.5512 | 降低 |
| normalized centroid distance | 0.3734 | 0.4372 | +0.0638 | 提升 |
| Fisher ratio | 0.1394 | 0.1912 | +0.0517 | 提升 |
| silhouette score | 0.0564 | 0.0916 | +0.0352 | 提升 |
| linear probe CV ROC-AUC | 0.8983 | 0.9068 | +0.0085 | 提升 |
| linear probe CV AP | 0.8232 | 0.8354 | +0.0122 | 提升 |

### MSH2 解读

MSH2 上，pLDDT-mask 后几乎所有 embedding 指标都朝预期方向变化：类间中心距离增大、类内方差降低、normalized centroid distance 与 Fisher ratio 上升，silhouette score 也明显提高。轻量线性 probe 的 CV ROC-AUC 和 AP 同样提升。

这说明在 MSH2 这个独立全长蛋白上，pLDDT-aware masking 不仅提高了 zero-shot scoring 表现，也使突变位点 embedding 更适合区分 pathogenic 与 benign。

## 与主实验的关系

该机制分析与前面实验形成互补：

| 实验 | 支撑的结论 |
|---|---|
| TP53 zero-shot | pLDDT-aware masking 在主蛋白上提升预测性能 |
| MSH2 zero-shot | pLDDT-aware masking 可跨蛋白全长复现 |
| Masking 严格消融 | 提升来自 pLDDT 选择位置，而不是随机 token 删除 |
| Embedding 空间分析 | masking 改变突变表征空间，使类别几何分离更清晰，尤其在 MSH2 上最一致 |

## 总体结论

Embedding 空间分析提供了一个机制层面的解释：pLDDT-aware masking 可能通过抑制低/中置信度结构 token 的噪声，使 SaProt 在突变位点形成更有判别性的表征。

具体而言：

| 观察 | 解释 |
|---|---|
| TP53 和 MSH2 上 centroid distance 均增大 | pathogenic 与 benign 的类别中心更远 |
| TP53 和 MSH2 上 Fisher ratio 均上升 | 类间分离相对于类内散布更强 |
| TP53 和 MSH2 上 silhouette score 均上升 | 两类在 embedding 空间中的标签一致性更好 |
| MSH2 上 within-class variance 降低，linear probe 提升 | masking 后表征空间更紧凑且更可线性区分 |
| TP53 上 linear probe 未提升 | 机制证据需保守，不能声称所有表征指标普遍提升 |

因此，本实验支持以下保守结论：

```text
pLDDT-aware masking 不仅改变 SaProt 的输出分数，也会改变突变位点的 embedding 空间。在 TP53 和 MSH2 上，masking 后 pathogenic 与 benign 的几何分离指标整体提高，尤其在 MSH2 上表现为类间距离增大、类内方差下降和线性 probe 表现提升。这为 masking 有效提供了机制层面的支持，但该证据应作为辅助解释，而不是替代 zero-shot bootstrap 结果的主证据。
```

## 推荐报告表述

```text
To further understand why pLDDT-aware masking improves zero-shot variant effect prediction, we analyzed the mutation-site embeddings extracted from SaProt before and after masking. For each variant, we masked the amino-acid token at the mutated position and extracted the last-layer hidden representation at that position. On both TP53 and MSH2, pLDDT-aware masking increased the pathogenic-benign centroid distance, normalized centroid distance, Fisher ratio, and silhouette score. On MSH2, masking also reduced within-class variance and improved cross-validated linear-probe AUC/AP. These results suggest that pLDDT-aware masking makes the mutation-site representation more discriminative by suppressing noisy low-confidence structure tokens while preserving reliable structural information. We treat this embedding analysis as supportive mechanistic evidence, with the primary performance claims still grounded in paired bootstrap zero-shot results.
```
