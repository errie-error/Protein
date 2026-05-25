# 🧬 《计算生物学》期末研究课题立项报告

## 📌 项目名称
**超越静态结构：基于结构置信度感知语言模型（SaProt）的突变效应预测与泛化性研究**
*(Beyond Static Structures: pLDDT-aware Structure-Conditioned Language Modeling for Mutation Effect Prediction)*

## 📖 1. 研究背景与核心动机 (Motivation)

在《计算生物学》课程中，我们学习了海量基因变异对表型的影响（第2讲），以及 AlphaFold2 在蛋白质结构预测上的革命性突破（第4、5讲）。然而，学术界最新共识指出：**AF2/ESMFold 等模型倾向于预测蛋白质稳定的“基态构象”（Ground-state conformation）**，并非为预测单点突变带来的微小自由能（$\Delta\Delta G$）变化而设计。因此，直接比较野生型（WT）与突变型（Mutant）的 3D 结构差异，信号往往过于微弱（即“突变盲区”）。

**本项目的核心创新在于：** 我们不强求 AI 去预测突变体的静态结构，而是引入最新的**结构感知蛋白质大模型（SaProt, ICLR 2024）**。我们将野生型的高精度 3D 结构作为环境约束（Condition），计算突变氨基酸在该 3D 几何环境下的**“结构条件突变对数似然比（Structure-conditioned mutational log-likelihood ratio, LLR）”**。此外，我们针对 AlphaFold 预测结构中存在的天然无序区（低 pLDDT 区域）导致的“结构噪音”问题，提出了基于置信度掩蔽的消融优化方案。

## 🎯 2. 研究对象与数据清洗 (Target Selection & Data Curation)

为保证模型能学习到真实的物理/生物学规律，而非死记硬背长序列中的无序噪音，本项目进行了极其严谨的研究对象收缩与数据清洗。

*   **研究对象：全长 TP53 肿瘤抑制蛋白**
    *   *弃用全长 BRCA1 的理由：* 全长超 1800 aa，存在大量低置信度无序区（IDR），提取结构 Token 易引入海量噪声，且容易爆显存。
    *   *选择全长 TP53 的理由：* 人类癌症中突变频率最高的基因，既包含高置信核心结构区，也包含低置信柔性区域，适合同时评估结构先验收益与 pLDDT 掩蔽策略。跨基因泛化测试将备选 `PTEN`。
*   **严谨的数据清洗 Pipeline (ClinVar)：**
    *   提取 ClinVar 中关于 TP53 的变异记录。
    *   **精准过滤：** 仅保留 `single amino-acid substitution`（单氨基酸错义突变）。
    *   **降噪：** 剔除所有 `VUS`（临床意义未明）、`Conflicting`（存在争议）的变异，仅合并明确的 `Pathogenic/Likely pathogenic` (正样本) 和 `Benign/Likely benign` (负样本)。
    *   预期清洗后获得极为纯净的数百至一千余条高质量二分类数据集。

## 🧪 3. 核心实验设计 (Experimental Design - The 3-Tier Architecture)

本项目将突破“调包微调”的局限，设计三层递进的硬核实验，全面拆解大模型的黑盒。

### 🔬 实验一：Zero-shot 基准评估（零样本结构条件 LLR）
*   **方法：** 不使用任何标签数据进行训练。直接将突变位点掩蔽（`[MASK]`），比较引入结构 Token 的 `SaProt (Seq + 3Di)` 与纯序列模型 `ESM-2 (Seq Only)` 的零样本打分（Zero-shot LLR）。
*   **目的：** 验证在无监督情况下，引入 3D 结构先验信息对致病突变识别能力的绝对提升（输出 ROC 曲线与 AUC 值对比）。

### 🔬 实验二：硬核消融实验（基于 pLDDT 的结构噪音掩蔽）
*   **假说：** 老师在课上强调了“静态照片”无法描述“动态视频”。我们认为，AlphaFold 给出的低置信度（pLDDT < 70）柔性 Loop 区域，其提取的 3Di 结构 Token 是不可靠的假阳性噪音，会误导语言模型。
*   **实验组别：**
    *   *Group A:* SaProt + 全序列完整 WT 3Di Token。
    *   *Group B (Ours):* SaProt + **Masked 3Di Token**（编写脚本，根据 AF2 的 pLDDT B-factor 列，将 pLDDT < 70 区域的结构 Token 替换为未知符号，仅保留高刚性区域的结构信息）。
*   **目的：** 证明剔除低置信度的动态区域结构噪音，能进一步提升模型的 LLR 预测准确率。这是本项目的**核心学术发现点**。

### 🔬 实验三：防泄露的监督微调（Spatial Splitting Fine-tuning）
*   **防弹级设计：** 传统的 8:2 随机划分极易导致“数据泄露（Data Leakage）”（模型记住了某个残基位点经常致病，而非学到了物理规律）。
*   **我们的方案：** 采用 **空间划分（Spatial Splitting）** 或 **跨基因划分**。使用 KMeans 基于 PDB 三维坐标将 TP53 突变分为左右两个半球，左半球训练，右半球测试；或在 TP53 上采用 LoRA（低秩微调）/ Linear Probing（线性探测），在 PTEN 上进行验证。
*   **目的：** 证明我们的模型真正学到了蛋白质空间不兼容的物理规律，具有极其强大的泛化鲁棒性。

## 👨‍💻 4. 四人团队分工矩阵 (Team Roles)

| 角色 | 姓名 | 核心职责与任务细节 |
| :--- | :--- | :--- |
| **数据科学家**<br>(Data Scientist) | [同学A] | 负责 ClinVar 数据库的清洗，严格执行单点错义突变的筛选逻辑。编写 Python 脚本解析 PDB 文件，提取 pLDDT 数值，并基于坐标完成 Spatial Split 空间数据集划分。 |
| **算法工程师**<br>(AI Engineer) | [同学B] | 负责 HuggingFace 环境部署。编写基于 SaProt 和 ESM-2 的 Zero-shot LLR 计算代码，生成突变评分。负责微调（LoRA/Linear Probing）的代码实现。 |
| **生信工程师**<br>(Bioinfo Engineer) | [同学C] | 负责干实验 Pipeline 搭建。调用 AlphaFold/ESMFold 获取野生型高精度 PDB，部署 `Foldseek`，完成 1D 氨基酸序列到双模态 `(Seq, 3Di)` 格式的自动化转换映射。 |
| **可视化与汇报**<br>(Visualizer & PM) | [同学D] | 负责实验数据汇总，使用 Python (Seaborn/Matplotlib) 绘制高标准 ROC 曲线及消融实验对比图。使用 **PyMOL** 在 3D 结构上高亮展示 pLDDT 掩蔽区域及模型判别结果，制作期末答辩 PPT。 |

## 🏆 5. 预期交付物 (Deliverables)

1.  **极具学术价值的发现：** 首次量化展示在突变预测中，“屏蔽 AlphaFold 低置信度区域的结构特征”能带来模型性能的提升。
2.  **高质量可视化：** PyMOL 渲染的 TP53 3D 突变热力图（结合空间划分数据集展示）。
3.  **开源级代码库：** 包含一键运行的 Data Cleaning、Foldseek Tokenization、Zero-shot LLR 和 LoRA Fine-tuning 的完整 Jupyter Notebook / Python 脚本。
4.  **期末研究报告：** 对标 Bioinformatics / ICLR 级别论文格式的详细学术报告。

## 📅 6. 项目进度安排 (Timeline)
*   **Week 1-2：** 环境配置（PyTorch, HuggingFace, Foldseek），完成 TP53 ClinVar 数据下载与严谨清洗。
*   **Week 3-4：** 获取 TP53 野生型 PDB，完成 3Di 词表转换；跑通 SaProt 与 ESM-2 的 Zero-shot LLR 并在全集上对比 AUC。
*   **Week 5-6：** 提取 pLDDT 分数，开展**核心消融实验（Masked 3Di vs Full 3Di）**，统计并分析结果。
*   **Week 7-8：** 划分空间验证集（Spatial Split），在顶层加入分类头进行轻量级微调（Linear Probing），进行性能测试。
*   **Week 9-10：** PyMOL 三维渲染可视化，整合图表，撰写期末报告与答辩 PPT，准备预演。

---
**附：核心技术栈与依赖库**
*   **计算框架：** `PyTorch`, `Transformers (HuggingFace)`, `PEFT (LoRA)`
*   **生信工具：** `Foldseek`, `Biopython`, `DSSP`, `PyMOL`
*   **基座模型：** `westlake-repl/SaProt_650M_AF2`, `facebook/esm2_t33_650M_UR50D`
