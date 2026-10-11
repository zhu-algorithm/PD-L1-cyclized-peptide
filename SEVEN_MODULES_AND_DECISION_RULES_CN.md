# 七模块决策规则与科学边界

本平台将“高亲和力、蛋白酶稳定、膜通透、口服暴露、体内靶点覆盖与可化学合成”拆解为可审阅的筛选字段。当前实现均是基于标准氨基酸序列组成的透明 **triage surrogate**，用于把候选送入下一步实验队列；并非训练完成的 PK/PD 模型，也不能替代实验数值。

| 模块 | 当前字段 | 如何进入决策 | 必需验证 |
|---|---|---|---|
| 1 靶点/文献 | PD-L1、PDB 5J89、His66/Gln116/Tyr123 | 定义受体背景与对接口袋 | 结构版本、文献证据与受体制备审查 |
| 2 环肽生成 | 标准 AA 序列、HELM、head-to-tail 概念环化 | 生成 5–20 aa 变体 | 合成路线和真实环化产物确认 |
| 3 亲和力/选择性 | `pd_l1_binding_proxy`、`affinity_triage_tier`、`selectivity_score` | `high-triage` 仅代表代理值 ≥0.75；不得写作高亲和力结论 | SPR/ITC/MST、PD-1/PD-L1 阻断实验、反筛 |
| 4 可开发性 | `protease_stability_score`、`permeability_score`、`oral_exposure_score`、毒性风险 | 联合筛除不利的序列组成；不直接给出半衰期、Papp 或 F% | 蛋白酶/血清稳定、Caco-2/MDCK、肠肝代谢和 PK |
| 5 分子对接 | `docking_score_kcal_mol` | 当前为确定性占位接口，不能作为 Vina/MD 结果 | 统一 PDB、网格、种子和 redocking 验证 |
| 6 候选排序 | `in_vivo_target_coverage_proxy`、可合成性、综合 priority | target coverage 代理组合结合、口服暴露、选择性和毒性；不是体内占有率 | 暴露量、游离浓度、组织分布、靶点占有/PD 标志物 |
| 7 实验闭环 | assay/value/unit 记录 | 保存为模型校准队列 | 批次、重复、质控与再训练方案 |

## 化学合成判断

`synthesizability` 是线性 SPPS 的序列复杂度代理。输出 `synthesis_risk_flags` 提示长肽、较高 Cys 比例、多 Trp 等可能需要工艺审查的因素。它不代表收率、纯度、成本或放大生产可行性；真实项目应进行 Fmoc SPPS 路线设计、粗品 LC-MS、纯化与环化收率评估。

## 推荐推进门槛

候选应同时满足：高优先级代理分数、无明显合成风险、可接受的稳定性/渗透/口服暴露代理，并完成已验证 docking QC 后，才进入合成与实验。最终“高亲和力”“口服可用”“体内靶点覆盖”只能由对应实验支持。
