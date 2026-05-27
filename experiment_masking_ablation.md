# 补充实验二：Masking 策略的严格消融

## 实验目的

本实验用于验证 pLDDT-aware SaProt masking 的提升是否真正来自“选对了需要掩蔽的位置”，而不是来自以下混杂因素：

| 可能混杂因素 | 需要排除的问题 |
|---|---|
| 减少了结构 token 数量 | 是不是随便 mask 一些 3Di token 也能提升？ |
| 改变了结构 token 分布 | 是不是只要改变输入形式就会提升？ |
| 保留低置信度结构信息 | 如果反过来保留低置信度结构 token、mask 高置信度 token，会不会也有效？ |

主实验使用 `pLDDT < 90` 作为 masking 阈值，因为前面的 TP53 和 MSH2 实验中，threshold 90 是最稳定、效果最强的设置。

## 对照组设计

| 策略 | 描述 | 作用 |
|---|---|---|
| SaProt full | 不进行 masking，使用完整 amino-acid + 3Di 序列 | 基线 |
| SaProt pLDDT-mask | mask 掉 pLDDT < 90 的结构 token | 实验组 |
| SaProt random-mask | 随机 mask 与 pLDDT-mask 相同数量的结构 token | 排除“只是减少 token 数量”的解释 |
| SaProt high-mask | 反向 mask pLDDT > 90 的高置信度结构 token，保留低/中置信度结构 token | 验证高置信度结构 token 是否关键 |

需要注意的是，`high-mask` 的解释应当保守：它主要说明“高置信度结构 token 很重要，不能被随意删除”，而不是说明“低置信度 token 完全没有信息”。

## 实现方式

在 `protein_project/structure.py` 中，SaProt 序列构建逻辑被扩展为同时生成三种 masking 序列：

| Payload key | 含义 |
|---|---|
| `masked_combined_seq` | pLDDT-aware masking，即 mask pLDDT < threshold 的结构 token |
| `random_masked_combined_seq` | 随机 masking，mask 数量与 pLDDT-aware masking 完全相同 |
| `high_masked_combined_seq` | 反向 masking，mask pLDDT > threshold 的高置信度结构 token |

在 `scripts/run_zero_shot.py` 中，zero-shot 评分脚本会自动检测并评分这些新增序列，输出以下分数列：

| 分数列 | 含义 |
|---|---|
| `saprot_full_score` | SaProt full |
| `saprot_masked_score` | pLDDT-aware masked SaProt |
| `saprot_random_masked_score` | random-mask SaProt |
| `saprot_high_masked_score` | high-mask SaProt |

在 `protein_project/benchmarks.py` 中，pairwise bootstrap comparison 增加了以下对比：

| Pairwise comparison | 回答的问题 |
|---|---|
| pLDDT-mask - random-mask | pLDDT 选点是否优于同数量随机 mask？ |
| pLDDT-mask - high-mask | 保留高置信度结构 token 是否重要？ |
| random-mask - full | 随机删除结构 token 是否本身有效？ |
| high-mask - full | 删除高置信度结构 token 是否会破坏性能？ |

运行命令如下：

```bash
PYTHONPATH=. python scripts/prepare_tp53_dataset.py --config configs/tp53_plddt90.yaml --foldseek third_party/foldseek/bin/foldseek
PYTHONPATH=. python scripts/run_zero_shot.py --config configs/tp53_plddt90.yaml

PYTHONPATH=. python scripts/prepare_tp53_dataset.py --config configs/msh2_plddt90.yaml --foldseek third_party/foldseek/bin/foldseek
PYTHONPATH=. python scripts/run_zero_shot.py --config configs/msh2_plddt90.yaml
```

随机 masking 使用固定 seed `0`，保证结果可复现。

## Mask 数量检查

| 蛋白 | pLDDT threshold | pLDDT-mask 数量 | random-mask 数量 | high-mask 数量 |
|---|---:|---:|---:|---:|
| TP53 | 90 | 186 | 186 | 206 |
| MSH2 | 90 | 503 | 503 | 428 |

其中 random-mask 与 pLDDT-mask 的 mask 数量完全一致，因此二者差异可以用于判断“mask 位置选择”是否重要。

## TP53 结果

TP53 数据集包含 373 个 ClinVar missense variants，其中 pathogenic/likely pathogenic 为 228 个，benign/likely benign 为 145 个。

### 总体指标

| Method | ROC-AUC | AUC 95% CI | AP | AP 95% CI |
|---|---:|---:|---:|---:|
| ESM-2 | 0.9250 | [0.8989, 0.9500] | 0.9565 | [0.9385, 0.9721] |
| SaProt full | 0.9623 | [0.9443, 0.9789] | 0.9776 | [0.9656, 0.9879] |
| SaProt pLDDT-mask | 0.9747 | [0.9602, 0.9876] | 0.9847 | [0.9751, 0.9926] |
| SaProt random-mask | 0.9411 | [0.9179, 0.9630] | 0.9672 | [0.9529, 0.9798] |
| SaProt high-mask | 0.4054 | [0.3456, 0.4703] | 0.5192 | [0.4687, 0.5802] |

### Pairwise bootstrap 差值

| Comparison | AUC diff | AUC diff 95% CI | AP diff | AP diff 95% CI |
|---|---:|---:|---:|---:|
| pLDDT-mask - full | +0.0124 | [+0.0017, +0.0244] | +0.0071 | [+0.0011, +0.0140] |
| pLDDT-mask - random-mask | +0.0336 | [+0.0147, +0.0543] | +0.0175 | [+0.0064, +0.0294] |
| pLDDT-mask - high-mask | +0.5693 | [+0.5066, +0.6318] | +0.4654 | [+0.4069, +0.5136] |
| random-mask - full | -0.0212 | [-0.0404, -0.0036] | -0.0104 | [-0.0211, -0.0002] |
| high-mask - full | -0.5569 | [-0.6169, -0.4959] | -0.4584 | [-0.5059, -0.4014] |

### TP53 结论

TP53 上，pLDDT-mask 显著优于 random-mask，AUC 差值为 `+0.0336`，95% CI 为 `[+0.0147, +0.0543]`；AP 差值为 `+0.0175`，95% CI 为 `[+0.0064, +0.0294]`。这说明 pLDDT-aware masking 的收益不是简单来自“mask 掉同样数量的结构 token”，而是来自“根据 pLDDT 选择了更合适的 masking 位置”。

同时，random-mask 显著低于 SaProt full，说明随机删除结构 token 不但不能解释提升，反而会破坏 SaProt 的结构信息。high-mask 的性能大幅下降，说明高置信度 3Di token 对 TP53 预测非常关键。

## MSH2 结果

MSH2 数据集包含 426 个 ClinVar missense variants，其中 pathogenic/likely pathogenic 为 141 个，benign/likely benign 为 285 个。

### 总体指标

| Method | ROC-AUC | AUC 95% CI | AP | AP 95% CI |
|---|---:|---:|---:|---:|
| ESM-2 | 0.8870 | [0.8482, 0.9216] | 0.7979 | [0.7219, 0.8665] |
| SaProt full | 0.9105 | [0.8747, 0.9421] | 0.8772 | [0.8193, 0.9228] |
| SaProt pLDDT-mask | 0.9278 | [0.8955, 0.9549] | 0.9125 | [0.8734, 0.9441] |
| SaProt random-mask | 0.8988 | [0.8610, 0.9304] | 0.8412 | [0.7703, 0.9008] |
| SaProt high-mask | 0.7104 | [0.6544, 0.7641] | 0.6255 | [0.5447, 0.7015] |

### Pairwise bootstrap 差值

| Comparison | AUC diff | AUC diff 95% CI | AP diff | AP diff 95% CI |
|---|---:|---:|---:|---:|
| pLDDT-mask - full | +0.0173 | [+0.0047, +0.0298] | +0.0354 | [+0.0068, +0.0715] |
| pLDDT-mask - random-mask | +0.0290 | [+0.0148, +0.0450] | +0.0714 | [+0.0302, +0.1191] |
| pLDDT-mask - high-mask | +0.2174 | [+0.1629, +0.2729] | +0.2870 | [+0.2155, +0.3573] |
| random-mask - full | -0.0117 | [-0.0273, +0.0026] | -0.0360 | [-0.0722, -0.0049] |
| high-mask - full | -0.2001 | [-0.2561, -0.1490] | -0.2517 | [-0.3175, -0.1908] |

### MSH2 结论

MSH2 上，pLDDT-mask 同样显著优于 random-mask。AUC 差值为 `+0.0290`，95% CI 为 `[+0.0148, +0.0450]`；AP 差值为 `+0.0714`，95% CI 为 `[+0.0302, +0.1191]`。这说明在独立的全长蛋白 MSH2 上，masking 位置选择仍然是关键因素。

high-mask 显著低于 SaProt full，说明删除高置信度结构 token 会明显损害性能。random-mask 相比 full 的 AUC 差值区间跨 0，但 AP 显著下降，因此随机 masking 不能作为有效替代策略。

## 总体结论

严格 masking 消融支持以下结论：

| 结论 | 证据 |
|---|---|
| pLDDT-aware masking 的收益不是来自简单 token 数量减少 | TP53 和 MSH2 上，pLDDT-mask 均显著优于同数量 random-mask |
| masking 位置选择很重要 | pLDDT-mask - random-mask 的 AUC/AP 差值在两个蛋白上均为正且 CI 不跨 0 |
| 高置信度结构 token 不能随意删除 | high-mask 在 TP53 和 MSH2 上均显著低于 SaProt full |
| 随机 masking 不是有效解释 | random-mask 在两个蛋白上均低于 pLDDT-mask，且在 TP53 上显著低于 full |

因此，本补充实验进一步支持项目核心假设：AlphaFold 低/中置信度区域的 3Di token 可能引入噪声，而基于 pLDDT 的选择性 masking 可以抑制这部分不可靠结构信息，同时保留高置信度结构 token 中的有用信号。

## 推荐报告表述

```text
为了排除性能提升仅来自减少结构 token 数量的可能性，我们设计了严格的 masking 消融实验。具体而言，我们将 pLDDT-aware masking 与同等 mask 数量的 random-mask 进行比较，并进一步加入反向的 high-mask 对照。结果显示，在 TP53 和 MSH2 两个蛋白上，pLDDT-aware masking 均显著优于 random-mask，说明提升来自对不可靠结构 token 的有选择性抑制，而不是简单的 token 删除。同时，mask 高置信度结构 token 的 high-mask 策略显著降低性能，表明高置信度 3Di token 仍然提供关键结构信息。因此，pLDDT-aware masking 的作用机制更合理地解释为：在保留可靠结构信息的同时，减少低/中置信度结构 token 带来的噪声。
```
