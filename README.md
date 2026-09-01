# PD-L1 环肽抑制剂发现平台

## LitRaPID-DT 接口

`POST /api/litrapid-export` 使用与 `/api/pipeline` 相同的参数，输出版本化的 `litrapid.pdl1-report.v1` 数据包。该数据包可直接提交给 [LitRaPID-DT](https://github.com/zhu-algorithm/LitRaPID) 的 `/api/pdl1-validation`，完成 RaPID 化学路线转换、多轮 mRNA display、RT-PCR/NGS 偏差模拟和最终候选排序。源头尾环化序列与 LitRaPID 生成的末端-Cys硫醚展示序列会分别保留，不能视为同一化学实体。

一个可运行的本地候选筛选界面，基于你提供的五阶段架构整理出可交互的“生成—评估—优先级排序”工作台。

## 页面导航

- `PDL1环肽筛选平台.html`：无需依赖、双击即可运行的 PD-L1 环肽筛选页面。
- `natural-products-ai-screening.html`：天然产物 AI 筛选平台的独立网页原型，涵盖谱图注释、去重复、来源追溯与新颖性优先级；不与环肽工作流混用。
- `source-code.html`：本地源代码导航页；上传 GitHub 后可直接通过仓库浏览这些文件。

## 八模块工作台

启动 `app.py` 后，浏览器首页为 `platform-dashboard.html`，提供以下可操作模块：

1. 靶点与文献知识库
2. 环肽序列生成
3. 结合亲和力与选择性预测
4. ADMET 与可合成性筛选
5. 分子对接接口
6. 多目标候选排序
7. 实验数据闭环学习记录
8. 抗体表位指导设计（Antibody Epitope-Guided Design）

对应的模块实现集中在 `platform_core.py`，网页接口在 `app.py`，界面在 `platform-dashboard.html`。

## 启动

在本目录执行：

```powershell
python app.py
```

浏览器会打开 `http://127.0.0.1:8765`。若未自动打开，请手动访问该地址。按 `Ctrl+C` 停止服务。

## 已实现

- HELM 与单字母环肽序列校验
- 基于种子序列的变体生成
- 可调多目标权重：拟合亲和力、渗透性、低毒性、药物相似性
- 候选排名、风险提示与 CSV 导出
- 每次运行的本地 JSON 留档（`runs/`）
- 基于可追溯抗体–PD-L1 复合物或用户提供接触残基列表的表位约束接口
- 表位兼容性代理分数，并将其纳入候选优先级排序

## 科学边界

为了保证可直接运行，本版不依赖 RDKit、PyTorch、对接工具或外部数据库；评分为透明的规则型 surrogate，仅适合候选的相对探索与流程演示。表位兼容性分数不等同于抗体表位重叠、竞争结合或 PD-1/PD-L1 阻断证据；生产环境必须导入有来源的抗体–PD-L1 结构或接触残基集，再经重对接和竞争实验验证。它不是经过训练和验证的 PD-L1 结合模型，也不能作为生物活性、ADMET 或临床决策依据。后续可将 `score_sequence` 替换为经过验证的 RDKit/ML/对接服务，并接入实验数据闭环。

