# Antibody Epitope-Guided Design（抗体表位指导设计）

## 目的

该模块将抗体–PD-L1 复合物的结构信息或经人工审查的接触残基列表，转化为环肽候选物的结构验证优先级。它用于提出“候选环肽可能覆盖或邻近阻断表位”的可检验假设。

## 输入要求

- 已注明来源的抗体–PD-L1 共晶或 cryo-EM 结构；或
- 与具体抗体、结构版本和残基编号绑定的人工审查接触残基列表。

默认的 `PDB:5N2F` 仅提供 PD-1/PD-L1 界面背景，不能代替特定抗体表位。生产使用时必须替换为指定抗体复合物的可追溯结构配置。

## 当前实现

`AntibodyEpitopeGuidedDesign` 在 `platform_core.py` 中输出：

- `epitope_reference`：表位配置、来源和接触残基说明；
- `antibody_epitope_compatibility_proxy`：基于序列组成与已配置表位的透明代理分数；
- `epitope_triage_tier`：是否优先进入表位重叠对接；
- `next_step`：结构对接和竞争实验的下一步要求。

候选排序新增 `epitope` 权重。该权重只用于决定结构验证队列，并不表示该候选物已经与任何抗体竞争结合。

## 内置的可切换结构配置

| 配置 ID | 靶点 | 设计区域 | 证据与用途 |
|---|---|---|---|
| `pdl1_antibody_pd1_facing_v_domain`（默认） | PD-L1 | N-terminal V-domain 的 PD-1-facing surface | PDB 5XXY（atezolizumab）和 5X8M（durvalumab）；用于 PD-L1 环肽候选排序。 |
| `pd1_antibody_loop_ensemble` | PD-1 | BC、CC′、C′D、FG loops 及邻近 PD-L1-facing surface | PDB 7WSL（dostarlimab）、7WVM（cemiplimab）与结构比较研究；仅用于平行的 PD-1 靶点设计，不能与 PD-L1 候选直接混排。 |
| `pdl1_patent_mapped_contacts` | PD-L1 | 报告的 QDAGVYRCMIS 区域及 D26/R113 上下文 | EP3455257B1、US10544225B2、WO2017055547A1 中披露的 Pepscan、HDX-MS、丙氨酸扫描与 epitope binning 证据；按专利逐项审查，不能视作普适表位。 |

调用 `run_pipeline(..., epitope_profile_id="pdl1_antibody_pd1_facing_v_domain")` 可选择配置；`GET /api/epitope/profiles` 会返回全部配置及来源链接。

## 验证路径

1. 使用与指定抗体一致的 PD-L1 构象进行重对接和姿势聚类；
2. 计算候选环肽和抗体的残基接触重叠、空间排斥或竞争可及性；
3. 通过竞争 SPR/BLI、PD-1/PD-L1 阻断实验或细胞功能实验进行验证；
4. 将已质控实验结果回写至实验数据闭环模块，更新后续排序。

## 科学边界

本模块不是抗体表位预测器，也不是竞争结合、亲和力、阻断活性或临床有效性的证据。所有当前数值均为本地规则型 surrogate，必须由可追溯结构、经验证的对接流程和独立实验支持。

