# TP53 ClinVar Mutation Effect Project

这个仓库实现了 `project.md` 对应的最小可运行版本，目标是围绕全长 `TP53` 做三件事：

- 清洗 `ClinVar` 的高置信 missense 变异
- 计算 `ESM-2` 与 `SaProt` 的 zero-shot 突变效应分数
- 在空间划分上做一个轻量级 `linear probe` 验证

## 目录结构

```text
configs/tp53.yaml
protein_project/
scripts/
project.md
requirements.txt
```

## 环境准备

```bash
python -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
```

如果你要跑 `SaProt`，还需要：

- 一个可用的 `Foldseek` 二进制
- 能访问 HuggingFace 模型 `westlake-repl/SaProt_650M_AF2`
- 能访问 HuggingFace 模型 `facebook/esm2_t33_650M_UR50D`

## 1. 构建数据

只做 `ClinVar + UniProt + AlphaFold` 数据构建：

```bash
python scripts/prepare_tp53_dataset.py --config configs/tp53.yaml
```

如果你同时已经安装好了 `Foldseek`，可以直接一并生成 `SaProt` 需要的结构感知序列：

```bash
python scripts/prepare_tp53_dataset.py --config configs/tp53.yaml --foldseek /path/to/foldseek
```

这个步骤会自动下载：

- `ClinVar variant_summary.txt.gz`
- `UniProt TP53 fasta`
- `AlphaFold TP53 pdb`

并在 `data/processed/` 下生成：

- `tp53_clinvar_clean.csv`
- `tp53_residue_table.csv`
- `tp53_saprot_sequences.json`（仅当传入 `--foldseek`）

## 2. 跑 zero-shot baseline 与 pLDDT ablation

```bash
python scripts/run_zero_shot.py --config configs/tp53.yaml
```

如果 `data/processed/tp53_saprot_sequences.json` 存在，这个脚本会同时输出：

- `esm2_score`
- `saprot_full_score`
- `saprot_masked_score`

输出文件：

- `results/tp53_zero_shot_scores.csv`
- `results/tp53_zero_shot_metrics.json`
- `results/tp53_zero_shot_roc.png`

## 3. 跑 spatial split linear probe

```bash
python scripts/run_linear_probe.py --config configs/tp53.yaml
```

输出文件：

- `results/tp53_linear_probe_metrics.json`
- `results/tp53_linear_probe_predictions.csv`

## 当前实现包含的关键逻辑

- `protein_project/data.py`
  - 下载原始数据
  - 解析 `ClinVar` 蛋白变异表示
  - 过滤到 `TP53`、missense、清晰 `Pathogenic/Benign` 标签
  - 与 `TP53` 参考序列对齐

- `protein_project/structure.py`
  - 从 `AlphaFold PDB` 提取 `CA` 坐标与 `pLDDT`
  - 调用 `Foldseek` 生成 `3Di`
  - 构造 `SaProt` 所需的 `combined sequence`
  - 基于三维坐标生成 `spatial_cluster`

- `protein_project/zero_shot.py`
  - `ESM-2` zero-shot LLR
  - `SaProt` zero-shot LLR
  - `pLDDT masked` 与 `full 3Di` 对比

- `protein_project/benchmarks.py`
  - `ROC-AUC` / `PR-AUC` / `Accuracy` / `F1`
  - `ROC` 绘图
  - `spatial split logistic regression`

## 你们答辩时可以直接讲的实验主线

1. `ESM-2` vs `SaProt full`
2. `SaProt full` vs `SaProt masked`
3. `Spatial split linear probe`

## 目前的边界

当前版本已经把主实验代码与数据流水线搭起来了，但还没有做这些扩展项：

- `PTEN` 跨基因泛化
- `LoRA` 微调
- `PyMOL` 自动作图
- `Bootstrap CI`
- 更细的 `ClinVar review status` 分层分析

如果后续需要，我建议按这个顺序继续补：

1. 先真的跑通 `prepare_tp53_dataset.py`
2. 再跑 `run_zero_shot.py`
3. 看结果后决定是否继续加 `PTEN` 或 `LoRA`
