# 研究方案

## 2023年全球结核病、2型糖尿病与矽肺负担重叠及共同宿主免疫机制：基于GBD 2023、WHO政策框架与跨疾病公共转录组的多层证据研究

**版本：V1.0**  
**日期：2026-08-07**  
**研究类型：全球横断面疾病负担分析 + 公共转录组二次分析 + 单细胞机制定位 + 基因调控网络 + 虚拟敲除**  
**GBD时间范围：仅使用2023年，不开展1990–2023趋势分析**

---

# 一、研究摘要

## 1.1 研究背景

结核病（tuberculosis, TB）仍是全球重要传染病。TB的发生、进展和治疗恢复受到宿主营养、免疫、代谢状态以及职业暴露共同影响。糖尿病和硅尘暴露/矽肺均已进入WHO结核病高风险因素与重点管理人群框架。WHO关于系统筛查的现行建议明确将硅尘暴露人群列为应系统筛查TB的重点人群，同时将糖尿病列为可结合当地流行病学特征优先考虑筛查的临床高风险状态。WHO 2025年TB与共病管理操作手册进一步将糖尿病纳入TB共病综合管理。

现有研究已经分别证明糖尿病与TB发生和治疗结局相关，硅尘暴露和矽肺增加TB风险，TB合并糖尿病存在明显宿主转录组异常，矽肺患者存在髓系细胞、巨噬细胞和干扰素相关免疫异常。公开数据库中尚未发现一个大规模人类转录组队列同时包含TB、T2DM和silicosis三种状态，因此本研究采用多层证据三角验证设计。

GBD 2023负责描述2023年全球三种疾病负担的空间重叠；WHO政策框架负责解释这些重叠热点的公共卫生意义；bulk RNA-seq和single-cell RNA-seq负责研究糖尿病和矽肺这两个不同TB风险状态是否在宿主免疫层面汇聚至共同的髓系细胞调控程序；基因调控网络和虚拟敲除进一步评估关键调控节点受到计算扰动后，共同异常程序是否出现逆转。

## 1.2 核心科学问题

### 问题1：全球分布

2023年全球哪些国家和地区同时具有较高的TB、T2DM和silicosis负担？

### 问题2：共同机制

糖尿病对TB宿主反应产生的附加转录扰动，与硅尘暴露背景下矽肺形成相关的转录扰动，是否存在稳定的共同基因、通路或共表达模块？

### 问题3：细胞定位

这些共同异常主要定位于哪些免疫细胞，尤其是否集中于单核细胞和巨噬细胞谱系？

### 问题4：计算扰动

对跨队列重复验证的核心调控基因进行虚拟敲除后，共同疾病相关调控网络和下游通路是否出现方向性逆转？

---

# 二、总体研究目标

## 2.1 总目标

构建一个由全球疾病负担、公共卫生政策、bulk transcriptomics、single-cell transcriptomics和gene regulatory network共同支持的跨疾病机制框架，阐明糖尿病和矽肺这两类TB高风险状态是否通过共同宿主免疫程序增加TB易感性、疾病进展或延迟治疗后的免疫恢复。

## 2.2 具体目标

1. 使用GBD 2023单一年份数据绘制全球TB、T2DM和silicosis负担地图。
2. 建立2023年三病负担重叠指数并识别triple-high-burden countries。
3. 将2023年高负担热点与WHO现行TB风险人群筛查和TB-DM共病管理框架对应。
4. 在GSE114192中定义糖尿病在TB背景下的附加转录效应。
5. 在GSE165489中定义硅尘暴露背景下矽肺形成相关的转录效应。
6. 在基因、通路、模块和调控网络层比较两个风险轴是否存在分子收敛。
7. 使用GSE181143、GSE193978、GSE193979和GSE249102进行TB-DM多队列验证。
8. 使用GSE264182进行矽肺肺局部bulk验证。
9. 使用GSE283452研究T2DM是否改变人巨噬细胞对M. tuberculosis刺激的应答。
10. 使用GSE174725、GSE326212和GSE192483完成单细胞和肺病灶定位。
11. 通过SCENIC/pySCENIC建立关键髓系细胞regulon。
12. 对高可信调控基因开展virtual knockout，并评估疾病signature reversal。

---

# 三、研究假设

## 3.1 人群层假设

H1：2023年全球TB、T2DM和silicosis负担存在显著空间异质性，并形成一组具有明显三病负担重叠的国家和地区。

H2：高silicosis负担和高T2DM负担同时存在的地区，其TB负担更可能处于高水平。

H3：三病高负担重叠地区与WHO关于silica exposure人群TB系统筛查及TB-DM综合管理的风险人群框架具有较强公共卫生相关性。

## 3.2 bulk transcriptomics层假设

H4：TB-DM与TB-only之间存在稳定的糖尿病附加宿主转录信号。

H5：silicosis与silica-exposed non-silicosis之间存在稳定的矽肺形成相关宿主转录信号。

H6：上述两种独立风险状态在若干免疫通路和调控模块上具有方向一致的收敛。

## 3.3 single-cell层假设

H7：共同异常信号主要富集于单核-巨噬细胞谱系。

H8：TB和silicosis肺部免疫微环境中存在可重复的macrophage/monocyte异常状态。

## 3.4 虚拟敲除层假设

H9：对跨bulk和single-cell重复验证的核心调控因子进行虚拟敲除，可使共同疾病相关表达程序向低风险或稳态方向移动。

---

# 四、总体研究结构

```text
GBD 2023，仅2023年
│
├── TB：年龄标化发病率
├── T2DM：年龄标化患病率
├── Silicosis：年龄标化患病率
└── 三病年龄标化DALY率敏感性分析
        │
        ▼
全球三病重叠负担与热点识别
        │
        ▼
WHO风险人群与筛查政策解释
        │
        ▼
Bulk RNA-seq主发现
│
├── GSE114192
│   TB-DM vs TB
│   → ΔDM|TB
│
└── GSE165489
    Silicosis vs silica-exposed non-silicosis
    → ΔSilicosis|Silica
        │
        ▼
跨疾病收敛分析
│
├── DEG direction
├── RRHO
├── GSEA
├── GO/KEGG/Reactome/Hallmark
├── WGCNA
└── PPI
        │
        ▼
多队列外部验证
        │
        ▼
Single-cell localization
│
├── GSE174725：Silicosis BALF
├── GSE326212：TB BAL
└── GSE192483：TB lung
        │
        ▼
Macrophage/monocyte regulon
        │
        ▼
SCENIC / GRN
        │
        ▼
Virtual knockout
        │
        ▼
跨疾病扰动方向一致性
```

---

# 五、GBD 2023全球横断面分析

## 5.1 数据来源

数据来源为Institute for Health Metrics and Evaluation的Global Burden of Disease Study 2023。

推荐从以下渠道提取：

- GBD Results Tool
- Global Health Data Exchange
- GBD 2023 cause/location codebook

本研究仅选择Year = 2023。

明确不开展：

- 1990–2022历史序列
- EAPC
- Joinpoint
- APC
- BAPC
- 年均变化百分比
- 时间趋势拐点
- 未来负担预测

GBD在本研究中的定位为2023年全球横断面风险景观。

## 5.2 疾病定义

### 5.2.1 TB

Primary cause：Tuberculosis，主分析使用all-form tuberculosis。

如GBD 2023允许进一步分层，可在补充材料描述drug-resistant TB或HIV-associated TB。主文保持TB总体定义。

### 5.2.2 T2DM

Primary：Type 2 diabetes mellitus。

选择T2DM的理由：

1. GEO主队列中的糖尿病主要为T2DM；
2. T2DM构成TB-DM共病的主要代谢背景；
3. 有利于GBD和转录组疾病定义保持一致。

Sensitivity：Diabetes mellitus total。

### 5.2.3 Silicosis

Primary：Silicosis。

Sensitivity：Pneumoconiosis total。

如GBD 2023提供下级cause，可补充描述coal workers' pneumoconiosis、asbestosis和other pneumoconiosis。主研究坚持silicosis，因为WHO对TB系统筛查的职业风险重点直接指向silica exposure。

## 5.3 主要GBD指标

### 5.3.1 TB主指标

\[
ASIR_{TB,2023}
\]

即2023年TB年龄标化发病率。

辅助指标：

- TB age-standardized DALY rate
- TB age-standardized mortality rate
- TB incidence number

### 5.3.2 T2DM主指标

\[
ASPR_{T2DM,2023}
\]

即2023年T2DM年龄标化患病率。

辅助指标：

- T2DM age-standardized DALY rate
- T2DM age-standardized incidence rate

### 5.3.3 Silicosis主指标

\[
ASPR_{Silicosis,2023}
\]

即2023年silicosis年龄标化患病率。

辅助指标：

- silicosis age-standardized DALY rate
- silicosis age-standardized mortality rate

## 5.4 三病指标标准化

TB使用发病率，T2DM和silicosis使用患病率，因此不直接相加绝对率。

### 方法A：百分位标准化

对每种疾病计算国家百分位数：

\[
P_{TB,i}=Percentile(ASIR_{TB,i})
\]

\[
P_{DM,i}=Percentile(ASPR_{DM,i})
\]

\[
P_{SIL,i}=Percentile(ASPR_{Silicosis,i})
\]

实际计算建议：

\[
P^*=\frac{rank-0.5}{N}
\]

避免0和1形成绝对边界。

### 共同负担指数

Primary：

\[
CBI_i=(P_{TB,i}\times P_{DM,i}\times P_{SIL,i})^{1/3}
\]

采用几何平均的原因是任一疾病负担明显偏低时，综合值会下降，更贴近三病重叠概念。

### 方法B：Z-score敏感性分析

三个疾病分别对log-transformed rate进行z标准化：

\[
Z_{combined}=\frac{Z_{TB}+Z_{DM}+Z_{SIL}}{3}
\]

## 5.5 triple-high-burden定义

Primary threshold：

- TB ≥ 全球第75百分位；
- T2DM ≥ 全球第75百分位；
- Silicosis ≥ 全球第75百分位。

三项同时满足定义为triple-high-burden country。

Sensitivity：

- 67th percentile
- 80th percentile
- 90th percentile

## 5.6 八类型国家分类

将三种疾病各自按high/low二分类，形成8类：

1. Low-Low-Low
2. TB high only
3. DM high only
4. Silicosis high only
5. TB + DM high
6. TB + Silicosis high
7. DM + Silicosis high
8. TB + DM + Silicosis high

该分类用于全球地图和政策解释。

## 5.7 无监督聚类

输入：

- standardized TB ASIR
- standardized T2DM ASPR
- standardized silicosis ASPR

优先使用PAM/k-medoids。

K选择：

- silhouette width
- gap statistic
- cluster stability

预期可能识别：

- high-TB dominant
- metabolic dominant
- occupational lung disease dominant
- triple-overlap
- low overall burden

聚类仅作描述性分类。

## 5.8 空间分析

分别计算：

- TB Moran's I
- T2DM Moran's I
- Silicosis Moran's I
- combined burden Moran's I

空间权重Primary：Queen contiguity。

Sensitivity：

- KNN k=4
- KNN k=6

局部热点：

- Local Moran's I
- Getis-Ord Gi*

输出high-high、low-low和spatial outlier。

## 5.9 横断面生态关联

该部分设为secondary/exploratory。

Outcome：

\[
Y_i=log(ASIR_{TB,i})
\]

Exposures：

\[
X_{1i}=log(ASPR_{T2DM,i})
\]

\[
X_{2i}=log(ASPR_{Silicosis,i})
\]

模型1：

\[
Y_i=\beta_0+\beta_1X_{1i}+\beta_2X_{2i}+\epsilon_i
\]

模型2：

\[
Y_i=\beta_0+\beta_1X_{1i}+\beta_2X_{2i}+\beta_3X_{1i}X_{2i}+\epsilon_i
\]

如获得GBD 2023 SDI：

\[
Y_i=\beta_0+\beta_1X_{DM}+\beta_2X_{Silicosis}+\beta_3X_{DM}X_{Silicosis}+\beta_4SDI+\epsilon_i
\]

建议增加GBD region或super-region fixed effects。

生态模型仅解释国家层面的负担共变，不作个体因果推断。

---

# 六、WHO政策框架

## 6.1 TB系统筛查

WHO现行TB systematic screening guidance将以下人群列为重点系统筛查群体：

- household and close contacts
- people living with HIV
- people exposed to silica
- people in prisons

糖尿病属于可结合当地流行病学和卫生系统条件优先考虑筛查的临床风险因素。

## 6.2 TB preventive treatment

WHO 2024 TB preventive treatment guideline覆盖包括silicosis和diabetes在内的多个高风险群体，为本研究的风险分层提供政策依据。

## 6.3 TB-diabetes共病管理

WHO 2025 Operational handbook on tuberculosis: module 6: tuberculosis and comorbidities纳入diabetes，强调早期识别、双向筛查、糖代谢管理和TB治疗过程中的综合照护。

## 6.4 本研究的政策分析定位

WHO部分采用policy relevance mapping。

GBD中识别的triple-high-burden countries用于讨论：

- TB筛查资源配置潜在重点
- occupational health与TB programme整合潜力
- diabetes care与TB screening协同价值

本方案不直接评价国家是否落实WHO建议。如后续获得国家级TB-DM、矿工TB筛查或职业健康政策实施指标，可扩展为policy implementation analysis。

---
# 七、bulk transcriptomics主分析

## 7.1 总体原则

所有数据集独立预处理、独立建模。跨平台、跨国家、跨疾病数据不直接拼成一个总表达矩阵。跨数据集整合在以下层面进行：

- effect size
- gene rank
- pathway
- co-expression module
- protein interaction network
- regulon
- cell state

## 7.2 主发现队列1：GSE114192

### 7.2.1 数据特征

- 疾病轴：TB-DM
- 组织：whole blood
- modality：bulk RNA-seq
- 样本量：249
- 组别：TB-only、TB-DM、TB-intermediate hyperglycaemia、DM-only、healthy
- 多地区来源

### 7.2.2 Primary contrast

\[
\Delta_{DM|TB}=Expression_{TB+DM}-Expression_{TB-only}
\]

该比较用于识别在TB背景下糖尿病额外引起的宿主表达变化。

### 7.2.3 Secondary contrasts

- TB-IH vs TB-only
- DM-only vs healthy
- TB-only vs healthy
- TB-DM vs DM-only

主线始终以TB-DM vs TB-only为核心。

## 7.3 GSE114192 metadata审计

必须建立sample-level metadata：

|变量|说明|
|---|---|
|sample_id|GSM或library ID|
|subject_id|受试者ID|
|group|TB/TBDM/TBIH/DM/HC|
|site|国家或研究中心|
|sex|如可用|
|age|如可用|
|HbA1c|如可用|
|BMI|如可用|
|treatment|治疗状态|
|batch|测序批次|
|library|library ID|

正式分析前检查：

1. 一个subject是否有多个library；
2. technical replicate；
3. missing group；
4. site与group是否完全混杂；
5. 是否存在治疗后样本；
6. 是否存在无法识别糖尿病定义的样本。

## 7.4 RNA-seq QC

如使用raw count：

- library size
- zero-count proportion
- expressed gene number
- PCA
- sample correlation
- sample-sample distance
- sex marker consistency
- outlier detection

低表达过滤优先使用edgeR `filterByExpr`，或采用CPM > 1且至少在最小组的一定比例样本中表达的规则。过滤标准在analysis_decisions.md中预先固定。

## 7.5 GSE114192统计模型

基础模型：

\[
Expression\sim Site+Group
\]

若metadata支持：

\[
Expression\sim Site+Sex+Age+Group
\]

如果多地区异质性明显，优先：

1. 每个site独立分析；
2. 提取gene-level beta和SE；
3. 后续random-effects meta-analysis。

---

# 八、主发现队列2：GSE165489

## 8.1 设计

- whole blood RNA-seq
- unexposed healthy：23
- silica-exposed non-silicosis：33
- silica-exposed silicosis：30

## 8.2 Primary contrast

\[
\Delta_{Silicosis|Silica}=Expression_{Silicosis}-Expression_{Silica\ exposed\ non-silicosis}
\]

该contrast共享硅尘暴露背景，主要用于识别矽肺形成或矽肺易感相关宿主变化。

## 8.3 Secondary contrasts

### 暴露效应

\[
\Delta_{SilicaExposure}=Expression_{Silica\ exposed\ non-silicosis}-Expression_{Unexposed}
\]

### 总疾病效应

\[
Expression_{Silicosis}-Expression_{Unexposed}
\]

两个secondary contrasts用于区分：

- silica exposure直接效应
- silicosis形成相关效应
- 从暴露阶段持续到疾病阶段的信号

---

# 九、差异表达分析

## 9.1 软件

RNA-seq：

- DESeq2
- edgeR
- limma-voom

Primary建议统一使用DESeq2或edgeR。

Microarray：limma。

## 9.2 保存字段

每个contrast至少保存：

- gene symbol
- Ensembl ID
- log2FC
- SE
- Wald/t statistic
- raw P
- FDR
- baseMean/average expression

## 9.3 DEG标准

Primary：

\[
FDR<0.05
\]

火山图和高可信候选筛选可加：

\[
|log_2FC|\ge0.5
\]

该效应阈值不用于ranked GSEA的基因过滤。

## 9.4 火山图

至少绘制：

1. GSE114192 TB-DM vs TB
2. GSE165489 silicosis vs exposed-no-silicosis
3. 关键外部验证队列

仅标记最终重复验证的hub genes、candidate regulators和少量top genes。

---

# 十、跨疾病基因层收敛分析

## 10.1 方向一致性

定义：

\[
logFC_{DM|TB}\times logFC_{Silicosis|Silica}>0
\]

包括两种情况：

- 两者均上调
- 两者均下调

## 10.2 候选分级

### Tier 1 shared genes

- 两个主数据集均FDR < 0.05
- effect direction一致
- 至少一个独立验证队列支持

### Tier 2 shared genes/modules

- gene-level显著性未同时达到FDR阈值
- effect direction稳定
- pathway/module层面重复
- single-cell支持

Venn diagram仅作为补充展示，主证据来自效应方向、rank和pathway concordance。

## 10.3 RRHO

使用Rank-Rank Hypergeometric Overlap比较两个全基因ranked lists。

推荐排序指标：

\[
RankScore=sign(logFC)\times -log_{10}(P)
\]

或直接使用Wald/t statistic。

关注：

- up-up overlap
- down-down overlap

同时展示up-down和down-up区域，用于发现机制分化。

---

# 十一、功能富集

## 11.1 GO

主分析采用GO Biological Process。

Cellular Component和Molecular Function放入补充分析。

同时进行：

- ORA
- ranked GSEA

## 11.2 KEGG

预设关注：

- tuberculosis
- phagosome
- lysosome
- NF-kappa B
- TNF
- HIF-1
- Toll-like receptor
- NOD-like receptor
- apoptosis
- glycolysis
- oxidative phosphorylation

结果排序以数据为准。

## 11.3 Reactome

重点观察：

- interferon signaling
- innate immune system
- neutrophil degranulation
- cytokine signaling
- complement
- extracellular matrix organization

## 11.4 Hallmark

建议将MSigDB Hallmark设为主要跨队列基因集体系之一，减少GO条目冗余。

预设关注：

- INTERFERON_GAMMA_RESPONSE
- INTERFERON_ALPHA_RESPONSE
- TNFA_SIGNALING_VIA_NFKB
- INFLAMMATORY_RESPONSE
- HYPOXIA
- GLYCOLYSIS
- OXIDATIVE_PHOSPHORYLATION
- COMPLEMENT
- APOPTOSIS

## 11.5 GSEA concordance

每个gene set分别获得：

\[
NES_{DM|TB}
\]

\[
NES_{Silicosis|Silica}
\]

方向一致定义：

\[
NES_1\times NES_2>0
\]

绘制NES scatter：

- 第一象限：共同激活
- 第三象限：共同抑制

并报告Spearman correlation和方向一致比例。

---

# 十二、WGCNA与模块保存

## 12.1 构网

GSE114192与GSE165489分别独立构建WGCNA网络。

不将两个表达矩阵直接合并。

## 12.2 模块关联

GSE114192：

\[
ModuleEigengene\sim TB-DM
\]

GSE165489：

\[
ModuleEigengene\sim Silicosis
\]

## 12.3 模块保存

使用modulePreservation评价：

- TB-DM相关模块是否在silicosis数据中保存
- silicosis相关模块是否在TB-DM数据中保存

常用解释：

- Zsummary < 2：保存证据弱
- 2–10：中等保存
- >10：较强保存

模块保存结果作为跨疾病网络层证据。

---

# 十三、PPI分析

## 13.1 输入基因

优先使用同时满足以下条件的候选：

- direction-concordant
- replicated
- pathway-supported

避免将单一数据集所有DEG直接输入STRING。

## 13.2 工具

- STRING
- Cytoscape
- cytoHubba
- MCODE

## 13.3 hub gene筛选

综合：

- degree
- MCC
- betweenness
- cross-dataset replication
- single-cell expression
- regulon evidence

最终核心hub genes控制在5–15个。

---

# 十四、TB-DM多队列外部验证

## 14.1 GSE181143

优势：

- 560 libraries
- India + Brazil
- baseline、month 2、month 6

主要用途：

1. TB-DM vs TB-only重复验证；
2. India和Brazil分别估计；
3. 检验site heterogeneity；
4. 评估共同module治疗后变化。

重复测量模型：

\[
Expression\sim Group\times Time+Site+(1|Subject)
\]

地区分层时：

\[
Expression\sim Group\times Time+(1|Subject)
\]

## 14.2 GSE193978

- 276 whole-blood RNA-seq
- baseline
- week 2
- month 2
- month 6
- TB-only
- TB-DM
- intermediate hyperglycaemia

主要模型：

\[
ModuleScore\sim Group\times Time+(1|Subject)
\]

核心问题：糖尿病相关共同免疫模块是否在抗TB治疗过程中消退更慢。

## 14.3 GSE193979

用于：

- good outcome vs poor outcome
- shared module与治疗结局关联
- DM是否改变shared module与结局之间的关系

如样本量支持：

\[
Outcome\sim SharedModule+DM+SharedModule\times DM+covariates
\]

## 14.4 GSE249102

组别：

- control
- TB
- TB-DM
- DM
- poorly controlled DM

用途：糖代谢控制梯度的方向性验证。

样本很小，仅进行：

- effect direction
- pathway/module score
- heatmap
- exploratory trend

不做复杂机器学习。

---
# 十五、GSE283452：T2DM改变巨噬细胞Mtb应答的机制桥梁

## 15.1 重要性

该数据集包含T2DM与non-T2DM人群的human alveolar macrophages和monocyte-derived macrophages，并包含live M. tuberculosis challenge和多个时间点。

它是本研究中连接T2DM和Mtb宿主反应最直接的人类实验型公开数据。

## 15.2 主要模型

\[
Expression\sim T2DM+Mtb+Time+CellType+T2DM\times Mtb+T2DM\times Time+\ldots
\]

核心项：

\[
T2DM\times Mtb
\]

代表糖尿病是否改变Mtb刺激引起的表达反应。

## 15.3 验证目标

检查前面发现的shared pathways：

- Mtb感染后是否被激活
- T2DM背景下是否增强或减弱
- alveolar macrophage与MDM是否一致
- 是否具有时间依赖性

重复时间点或同一供体不同细胞条件必须以donor-aware模型处理。

---

# 十六、矽肺肺局部bulk验证：GSE264182

组别包括：

- chronic cough control
- silica-exposed no disease
- simple silicosis
- complicated silicosis

构建有序疾病严重程度：

\[
Control\rightarrow Exposure\rightarrow Simple\rightarrow Complicated
\]

主要分析：

- trend test
- ordinal contrast
- shared module score
- ranked GSEA

样本量有限，优先报告效应量和方向一致性。

---

# 十七、single-cell分析总体原则

## 17.1 核心数据集

### Silicosis

GSE174725：

- BALF
- silica-exposed without silicosis
- silicosis
- 5 donors

### TB airway

GSE326212：

- BAL
- recent household contacts
- controllers
- progressors
- active TB

### TB lung

GSE192483：

- lung resection
- TB lesion
- adjacent less-involved tissue
- 6 patients，11 samples

### T2DM补充

GSE268210：

- PBMC
- T2DM vs healthy
- scRNA/scTCR/scBCR

## 17.2 single-cell QC

每个数据集根据自身分布确定QC阈值，不机械固定统一cutoff。

流程：

1. 查看nFeature分布；
2. 查看nCount分布；
3. 查看percent.mt；
4. 每个sample独立识别极端值；
5. doublet detection；
6. ambient RNA评估；
7. 识别低质量cluster；
8. 保存所有QC前后细胞数。

工具：

- Seurat或Scanpy
- DoubletFinder或scDblFinder
- SoupX或CellBender

## 17.3 批次整合

Harmony、CCA或scVI可用于：

- visualization
- clustering
- annotation

差异表达统计不直接依赖整合后的校正表达值。

Primary inferential strategy：donor-level pseudobulk。

---

# 十八、细胞注释

## 18.1 一级细胞类型

- macrophage
- monocyte
- neutrophil
- dendritic cell
- CD4 T
- CD8 T
- NK
- B
- plasma cell
- epithelial
- mast cell

## 18.2 巨噬细胞亚型

根据数据和已发表marker进一步区分：

- alveolar macrophage
- inflammatory macrophage
- monocyte-derived macrophage
- interferon-responsive macrophage
- lipid-associated macrophage
- fibrotic macrophage
- proliferating macrophage

亚型命名必须有marker、参考文献和cluster表达证据支持。

---

# 十九、single-cell统计分析

## 19.1 细胞组成

每个donor计算：

\[
Proportion_{celltype}=\frac{n_{celltype}}{n_{allcells}}
\]

候选方法：

- propeller
- scCODA
- Dirichlet regression

## 19.2 donor-level pseudobulk

聚合单位：

\[
Donor\times CellType
\]

Silicosis macrophage primary contrast：

\[
Silicosis\ vs\ ExposedNoSilicosis
\]

TB macrophage contrasts根据metadata预设：

- ActiveTB vs Controller
- Progressor vs Controller

## 19.3 shared module score

将bulk发现的shared genes构造成module。

每个cell计算：

- UCell
- AUCell
- AddModuleScore作为补充

Primary优先UCell或AUCell。

统计比较以donor为单位，不以cell为独立样本。

---

# 二十、共同细胞状态判定规则

一个高可信shared myeloid state至少满足：

1. TB-DM bulk支持；
2. silicosis bulk支持；
3. silicosis scRNA中macrophage/monocyte支持；
4. TB scRNA相同细胞谱系支持；
5. effect direction一致；
6. 至少一个独立bulk或肺组织验证队列支持。

结果若指向neutrophil、T cell或dendritic cell，则按数据调整主细胞类型。

---

# 二十一、SCENIC / pySCENIC

## 21.1 目的

从共同差异程序中识别：

- upstream transcription factor
- regulon
- TF-target relationship

## 21.2 分析对象

Primary：

- macrophage
- monocyte

每个主要细胞谱系独立分析。

## 21.3 流程

1. expression matrix准备；
2. GRN inference；
3. motif enrichment；
4. regulon construction；
5. AUCell regulon activity；
6. disease vs reference comparison；
7. cross-disease regulon concordance。

## 21.4 跨疾病regulon

对每个regulon计算：

\[
Effect_{Silicosis}
\]

和：

\[
Effect_{TB}
\]

方向一致定义：

\[
Effect_{Silicosis}\times Effect_{TB}>0
\]

再结合TB-DM bulk regulator evidence进行优先级排序。

---

# 二十二、虚拟敲除

## 22.1 候选基因进入标准

候选TF/regulator至少满足以下4项：

1. TB-DM bulk显著或稳定方向；
2. silicosis bulk显著或稳定方向；
3. TB-DM外部队列重复；
4. silicosis/TB single-cell中定位到相同细胞类型；
5. SCENIC regulon显著；
6. PPI/GRN中心性较高。

严格候选数量建议5–10个。

## 22.2 工具

Primary：scTenifoldKnk。

备选：

- CellOracle
- 其他gene regulatory network perturbation framework

## 22.3 构网策略

对关键细胞类型分别建立：

- disease GRN
- reference GRN

降低细胞数量造成伪精确的措施：

- 每个donor等量抽样
- bootstrap
- leave-one-donor-out
- 多次重复构网
- 提取稳定edge

## 22.4 KO输出

每个KO保存：

- downstream perturbed genes
- perturbation magnitude
- pathway enrichment
- shared-module response

## 22.5 disease-signature reversal

定义shared disease signature：

\[
S=\{g_1,\ldots,g_k\}
\]

保留每个基因疾病方向。

虚拟KO产生perturbation vector：

\[
K=\{\Delta g_1,\ldots,\Delta g_k\}
\]

计算：

\[
\rho=Spearman(S,K)
\]

\[
\rho<0
\]

表示KO扰动与疾病signature整体方向相反。

再进行GSEA：

- disease-up genes在KO后是否获得negative NES
- disease-down genes在KO后是否获得positive NES

## 22.6 跨疾病KO验证

最高优先级candidate regulator需满足：

- silicosis macrophage网络KO后产生逆转；
- TB macrophage网络KO后产生相同方向逆转；
- downstream pathway具有一致性。

虚拟KO输出属于计算预测，不用于宣称实验因果或已确认治疗靶点。

---

# 二十三、细胞通讯分析

该部分为可选secondary analysis。

工具：

- CellChat
- CellPhoneDB
- LIANA

仅在以下条件满足后执行：

- donor数和组间结构足以比较；
- shared cell state明确；
- shared module已经稳定重复。

重点可观察：

- macrophage-T cell
- macrophage-epithelial
- monocyte-macrophage
- TNF
- IFN
- CXCL
- IL1
- TGF-beta signaling

---

# 二十四、跨数据集Meta-analysis

## 24.1 TB-DM gene-level meta

各TB-DM队列独立估计：

\[
\hat\beta_j,SE_j
\]

使用random-effects meta-analysis。

输出：

- pooled effect
- 95% CI
- tau-squared
- I-squared

## 24.2 异质性处理

优先分队列估计后meta，因为不同数据集存在：

- 国家差异
- ancestry差异
- library差异
- sequencing platform差异
- clinical definition差异
- treatment phase差异

如I-squared较高，进一步按国家、治疗时间或糖尿病控制状态分层。

---

# 二十五、候选基因最终分级

## Tier A

满足：

- 两个主bulk风险轴方向一致；
- TB-DM至少2个外部队列验证；
- single-cell定位一致；
- regulon或PPI支持；
- virtual KO出现疾病signature逆转。

Tier A构成最终核心调控候选。

## Tier B

满足：

- bulk方向一致；
- 至少一个验证队列；
- single-cell支持。

## Tier C

仅在单一数据集显著或只有网络中心性证据，无稳定外部重复。

Tier C不进入主要机制结论。

---

# 二十六、阴性结果决策规则

## 26.1 共同DEG数量很少

继续进行：

- RRHO
- ranked GSEA
- module preservation
- regulon concordance

跨疾病机制评价不依赖完全相同的单基因列表。

## 26.2 两个风险轴方向不同

若大部分核心通路呈相反方向，研究转向divergent host programs，分析两类TB风险状态通过不同免疫路径影响TB的可能性。

## 26.3 shared signal未定位到macrophage

根据真实结果转向neutrophil、T cell、dendritic cell或其他主要细胞谱系。

## 26.4 virtual KO无明显逆转

报告共同signature可能由多节点共同维持，不继续无限筛选候选基因以追求阳性结果。

---

# 二十七、敏感性分析

## 27.1 GBD

1. triple-high threshold：75th → 67th/80th/90th percentile；
2. percentile composite → z-score；
3. silicosis → pneumoconiosis total；
4. T2DM → diabetes total；
5. disease-specific primary rates → unified age-standardized DALY rates；
6. Queen adjacency → KNN spatial weights；
7. inclusion/exclusion of small island states；
8. combined index geometric mean → arithmetic mean。

## 27.2 Bulk

1. DESeq2 → edgeR；
2. different low-expression filters；
3. site-stratified analysis；
4. remove influential outliers；
5. leave-one-site-out；
6. effect-size threshold sensitivity。

## 27.3 Single-cell

1. UCell vs AUCell；
2. pseudobulk model alternatives；
3. leave-one-donor-out；
4. different annotation granularity；
5. rare-cell exclusion；
6. integration method sensitivity。

## 27.4 GRN

1. different random seeds；
2. donor-balanced subsampling；
3. SCENIC vs alternative GRN；
4. candidate KO rank stability；
5. bootstrap consensus network。

---

# 二十八、多重比较与统计显著性

- DEG：Benjamini-Hochberg FDR
- GSEA：FDR
- GO/KEGG/Reactome：FDR
- regulon：FDR
- cell-type comparisons：FDR
- spatial local tests：适当FDR校正

主要结果同时报告效应量和95% CI，避免仅依赖P值。

---

# 二十九、缺失数据与数据审计

GBD：

- 仅使用GBD正式estimate；
- 不自行插补国家疾病负担。

GEO：

- group不可识别样本排除对应primary contrast；
- demographic covariate缺失不做简单均值插补；
- 主模型优先使用完整且跨组可比的协变量；
- 每个排除样本记录原因。

---

# 三十、Primary与Secondary outcomes

## 30.1 Population primary outcomes

- 2023 triple-high-burden classification
- combined burden percentile index

## 30.2 Transcriptomic primary outcomes

- direction-concordant gene program
- concordant Hallmark/Reactome pathways
- replicated co-expression modules

## 30.3 Single-cell primary outcome

共同bulk module是否在TB和silicosis中定位至相同髓系细胞谱系。

## 30.4 Perturbation primary outcome

virtual KO后shared disease signature是否出现方向性逆转。

---
# 三十一、预设关键生物学方向

以下通路仅作为prior biological interest，不作为强制阳性结果：

1. IFN-gamma signaling
2. TNF/NF-kappa B
3. NLRP3 inflammasome
4. phagosome
5. lysosome
6. macrophage activation
7. neutrophil degranulation
8. complement
9. oxidative stress
10. glycolysis
11. HIF-1 signaling
12. mitochondrial metabolism
13. extracellular matrix remodeling
14. TGF-beta
15. apoptosis/autophagy

---

# 三十二、统计软件与环境

## 32.1 R

建议R >= 4.4。

主要包：

- data.table
- tidyverse
- DESeq2
- edgeR
- limma
- fgsea
- clusterProfiler
- ReactomePA
- msigdbr
- WGCNA
- metafor
- ComplexHeatmap
- ggplot2
- sf
- spdep
- tmap
- Seurat
- muscat
- UCell
- AUCell
- CellChat

## 32.2 Python

建议Python >= 3.10。

主要包：

- scanpy
- anndata
- scvi-tools
- decoupler
- gseapy
- pySCENIC

## 32.3 Cytoscape

- STRING app
- cytoHubba
- MCODE

所有package version写入sessionInfo或environment.yml。

---

# 三十三、可复现项目目录

```text
TB_DM_SILICOSIS_PROJECT/
│
├── 00_protocol/
│   ├── protocol_v1.md
│   ├── analysis_decisions.md
│   └── deviations_from_protocol.md
│
├── 01_GBD2023/
│   ├── raw/
│   ├── metadata/
│   ├── processed/
│   ├── scripts/
│   └── figures/
│
├── 02_GEO_bulk/
│   ├── GSE114192/
│   ├── GSE165489/
│   ├── GSE181143/
│   ├── GSE193978/
│   ├── GSE193979/
│   ├── GSE249102/
│   ├── GSE283452/
│   └── GSE264182/
│
├── 03_scRNA/
│   ├── GSE174725/
│   ├── GSE326212/
│   ├── GSE192483/
│   └── GSE268210/
│
├── 04_cross_disease/
│   ├── RRHO/
│   ├── GSEA/
│   ├── WGCNA/
│   ├── PPI/
│   └── meta/
│
├── 05_GRN_KO/
│   ├── scenic/
│   ├── scTenifoldKnk/
│   └── signature_reversal/
│
├── 06_results/
│   ├── tables/
│   ├── figures/
│   └── supplementary/
│
└── 07_manuscript/
```

---

# 三十四、推荐脚本执行顺序

```text
00_environment.R
01_download_GBD2023.R
02_clean_GBD2023.R
03_GBD2023_descriptive.R
04_GBD2023_hotspot.R
05_GBD2023_spatial.R
06_GBD2023_ecological.R

10_GEO_metadata_audit.R
11_GSE114192_QC.R
12_GSE114192_DE.R
13_GSE165489_QC.R
14_GSE165489_DE.R

20_cross_disease_RRHO.R
21_cross_disease_GSEA.R
22_cross_disease_WGCNA.R
23_cross_disease_PPI.R

30_TBDM_validation.R
31_silicosis_validation.R
32_longitudinal_validation.R
33_GSE283452_Mtb_T2DM.R

40_scRNA_GSE174725.R
41_scRNA_GSE326212.R
42_scRNA_GSE192483.R
43_scRNA_cross_disease.R

50_SCENIC.py
51_regulon_validation.R
52_virtual_knockout.R
53_signature_reversal.R

60_final_tables.R
61_final_figures.R
62_reproducibility_audit.R
```

---

# 三十五、预期主图

## Figure 1：Study design

展示GBD 2023、WHO、bulk discovery、external validation、scRNA、GRN和virtual KO之间的证据链。

## Figure 2：2023全球三病负担

A. TB ASIR map  
B. T2DM ASPR map  
C. Silicosis ASPR map  
D. triple-burden classification map  
E. combined burden index  
F. unsupervised cluster

## Figure 3：两个主要bulk风险轴

A. GSE114192 volcano  
B. GSE165489 volcano  
C. RRHO  
D. logFC concordance  
E. shared-gene heatmap

## Figure 4：Pathway convergence

A. Hallmark NES scatter  
B. Reactome/GO dotplot  
C. shared pathway network  
D. WGCNA module preservation

## Figure 5：外部验证

A. TB-DM gene/module forest plot  
B. GSE181143 India/Brazil  
C. GSE193978 longitudinal module trajectory  
D. GSE193979 outcome association  
E. GSE283452 T2DM × Mtb response

## Figure 6：Single-cell localization

A. GSE174725 UMAP  
B. GSE326212 UMAP  
C. shared module score  
D. macrophage pseudobulk  
E. shared myeloid state  
F. GSE192483 lung lesion validation

## Figure 7：Regulatory network and virtual KO

A. SCENIC regulon  
B. shared regulator network  
C. virtual KO perturbation  
D. disease-signature reversal  
E. TB vs silicosis KO concordance

## Figure 8：Integrated mechanism model

呈现silica exposure/silicosis和T2DM两条风险路径汇聚至共同或部分重叠的host immune program，并连接TB susceptibility/progression和WHO风险人群管理框架。

---

# 三十六、预期补充图

- GBD indicator distributions
- percentile threshold sensitivity
- Moran's I
- alternative spatial weights
- PCA and sample QC
- site-stratified volcano plots
- complete GSEA results
- WGCNA diagnostics
- full PPI networks
- scRNA QC
- cell marker heatmap
- donor-level pseudobulk sensitivity
- leave-one-donor-out results
- virtual KO bootstrap stability

---

# 三十七、预期主表

## Table 1

2023年TB、T2DM和silicosis负担最高国家与地区。

## Table 2

triple-high-burden countries及对应WHO风险人群政策意义。

## Table 3

所有纳入公共转录组数据集的样本、组织、平台、contrast和分析角色。

## Table 4

跨疾病重复验证的shared pathways/modules。

## Table 5

最终Tier A/Tier B candidate regulators。

---

# 三十八、创新点

## 创新1：职业肺病、代谢病和感染性疾病的统一宿主易感性框架

研究将silica/silicosis、T2DM和TB放在同一公共卫生和宿主免疫框架下评估。

## 创新2：GBD 2023仅承担当前全球空间定位

研究不复制常见的1990–2023趋势分析。GBD部分集中回答2023年全球风险重叠在哪里，为后续政策和机制分析提供人群背景。

## 创新3：条件差异设计

两个核心contrast为：

\[
TB-DM-TB
\]

和：

\[
Silicosis-SilicaExposedNoSilicosis
\]

前者控制TB疾病背景，后者控制silica exposure背景，优于三个疾病分别对healthy control后直接取交集。

## 创新4：跨疾病transcriptomic triangulation

通过effect direction、RRHO、GSEA、module preservation和regulon进行多层收敛验证。

## 创新5：肺部single-cell定位

silicosis BALF与TB BAL/BALF在组织环境上具有较高可比性，可用于确定共同异常属于哪类肺部免疫细胞。

## 创新6：virtual knockout形成可检验的机制预测

只有经过bulk、外部队列和single-cell重复的调控节点才进入虚拟敲除，从而形成发现、定位、调控和扰动的闭环。

---

# 三十九、主要局限

1. GBD分析属于国家层面横断面生态分析；
2. GBD不能识别同一患者是否同时患TB、T2DM和silicosis；
3. GBD估计受到各国原始数据质量和模型假设影响；
4. TB-DM和silicosis公共数据来自不同国家和人群；
5. whole-blood bulk expression受到细胞组成差异影响；
6. silicosis single-cell数据生物学重复数有限；
7. virtual KO属于计算扰动预测；
8. WHO policy mapping主要评价公共卫生相关性，缺少统一的国家执行程度指标；
9. 当前无单一公开队列同时包含TB、T2DM和silicosis；
10. 最终共享机制需要未来实验研究验证。

---

# 四十、结论语言边界

可以使用：

- cross-disease convergent host-response program
- shared immune pathway
- common macrophage-associated regulatory program
- overlapping population burden
- policy-relevant high-burden settings
- computational perturbation predicts

避免使用：

- triple-comorbidity molecular signature
- virtual knockout proves causality
- GBD proves diabetes causes TB
- GBD proves silicosis causes TB
- WHO policy implementation failure
- therapeutic target confirmed

---

# 四十一、论文主线

2023年全球TB、T2DM和silicosis负担具有明显空间异质性，部分国家呈现三病负担重叠。WHO现行TB政策分别将silica exposure和diabetes纳入风险管理框架。跨疾病公共转录组进一步检验糖尿病对TB宿主反应造成的附加扰动与矽肺形成相关宿主扰动是否在免疫通路和调控模块上发生收敛。单细胞分析定位这些共同程序的主要细胞来源。若结果支持预设假设，进一步聚焦macrophage/monocyte谱系。最终通过基因调控网络和虚拟敲除评估核心调控节点受到扰动后共同异常程序的可逆性。

---

# 四十二、候选论文题目

## 题目A

**Convergent host immune programs linking diabetes and silicosis to tuberculosis: evidence from GBD 2023, WHO risk-group policies, and cross-disease transcriptomics**

## 题目B

**Global overlap and shared host mechanisms of tuberculosis, type 2 diabetes, and silicosis: a GBD 2023 and cross-disease transcriptomic study**

## 题目C

**From global co-burden to macrophage regulatory programs: integrating GBD 2023 and human transcriptomics across tuberculosis, diabetes, and silicosis**

---

# 四十三、最低可发表版本

若single-cell或virtual KO因数据结构限制无法完成，最低完整版本保留：

1. GBD 2023三病全球重叠；
2. WHO政策框架；
3. GSE114192；
4. GSE165489；
5. RRHO；
6. GO/KEGG/Reactome/Hallmark GSEA；
7. GSE181143和GSE193978验证；
8. GSE174725或GSE326212至少一个single-cell定位。

---

# 四十四、推荐完整版本

1. GBD 2023横断面全球负担；
2. triple-high hotspot；
3. WHO policy relevance；
4. GSE114192 TB-DM primary contrast；
5. GSE165489 silicosis primary contrast；
6. RRHO；
7. GSEA；
8. WGCNA module preservation；
9. PPI；
10. GSE181143跨地区验证；
11. GSE193978 longitudinal recovery；
12. GSE193979 outcome；
13. GSE283452 macrophage Mtb challenge；
14. GSE174725 silicosis scRNA；
15. GSE326212 TB scRNA；
16. GSE192483 TB lung lesion；
17. SCENIC；
18. virtual KO；
19. KO signature reversal。

---

# 四十五、分析前审计清单

## GBD

- [ ] GBD Results Tool确认使用GBD 2023
- [ ] Year仅选择2023
- [ ] location层级统一
- [ ] sex主分析选择both sexes
- [ ] age选择age-standardized
- [ ] TB选择incidence rate
- [ ] T2DM选择prevalence rate
- [ ] silicosis选择prevalence rate
- [ ] 三病DALY sensitivity下载完成
- [ ] 95% uncertainty interval完整保留
- [ ] cause/location codebook保存

## Bulk

- [ ] 每个GEO逐GSM审计
- [ ] subject ID确认
- [ ] technical replicate确认
- [ ] site确认
- [ ] raw count与processed matrix来源记录
- [ ] gene ID转换版本记录
- [ ] low-expression规则固定
- [ ] outlier处理写入deviation log

## Single-cell

- [ ] donor ID确认
- [ ] library与donor关系确认
- [ ] cell数量不当作统计n
- [ ] doublet处理
- [ ] ambient RNA评估
- [ ] annotation marker保存
- [ ] pseudobulk设计矩阵检查

## GRN

- [ ] 候选基因来源可追溯
- [ ] 不从KO结果反向选择候选
- [ ] random seed保存
- [ ] donor-balanced subsampling
- [ ] bootstrap次数预先固定
- [ ] leave-one-donor-out完成

---

# 四十六、建议预注册的Primary analyses

建议在正式运行前预注册以下6项：

1. GBD triple-high threshold = 75th percentile；
2. GBD primary disease metrics = TB ASIR + T2DM ASPR + silicosis ASPR；
3. GSE114192 primary contrast = TB-DM vs TB-only；
4. GSE165489 primary contrast = silicosis vs silica-exposed non-silicosis；
5. pathway primary framework = Hallmark GSEA + Reactome GSEA；
6. single-cell primary inferential method = donor-level pseudobulk。

GO、KEGG、PPI、WGCNA、CellChat和virtual KO设为secondary/mechanistic analyses。

---

# 四十七、数据与政策来源

## GBD

Institute for Health Metrics and Evaluation. Global Burden of Disease Study 2023.

GBD Results Tool:  
https://vizhub.healthdata.org/gbd-results/

IHME:  
https://www.healthdata.org/

GHDx:  
https://ghdx.healthdata.org/

## WHO

WHO. Tuberculosis: Systematic screening.  
https://www.who.int/news-room/questions-and-answers/item/systematic-screening-for-tb

WHO consolidated guidelines on tuberculosis. Module 2: screening, systematic screening for tuberculosis disease.  
https://www.who.int/publications-detail-redirect/9789240022676

WHO consolidated guidelines on tuberculosis Module 1: prevention, tuberculosis preventive treatment, second edition.  
https://www.who.int/publications/i/item/9789240096196

WHO operational handbook on tuberculosis: module 6: tuberculosis and comorbidities, 3rd ed.  
https://www.who.int/publications/i/item/9789240103276

## GEO核心数据

GSE114192  
https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE114192

GSE181143  
https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE181143

GSE193978  
https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE193978

GSE193979  
https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE193979

GSE249102  
https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE249102

GSE283452  
https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE283452

GSE165489  
https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE165489

GSE264182  
https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE264182

GSE174725  
https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE174725

GSE326212  
https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE326212

GSE328391  
https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE328391

GSE192483  
https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE192483

GSE268210  
https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE268210

GSE184050  
https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE184050

---

# 四十八、最终执行建议

整篇研究保持以下单一逻辑链：

\[
2023\ Global\ overlap
\rightarrow
WHO\ risk-group\ relevance
\rightarrow
two\ independent\ host-risk\ contrasts
\rightarrow
cross-disease\ molecular\ convergence
\rightarrow
cellular\ localization
\rightarrow
regulatory\ perturbation
\]

篇幅建议：

- GBD 2023：约20%–25%
- transcriptomics + single-cell：约60%–65%
- WHO与公共卫生转化：约10%–15%

GBD负责回答当前全球负担在哪里。转录组负责回答两个TB风险状态共享什么宿主程序。single-cell负责回答该程序发生在哪类细胞。GRN和virtual KO负责回答哪些调控节点可能维持该程序以及计算扰动后会发生什么变化。

---

# 四十九、一句话研究问题

**在2023年全球疾病负担格局中，TB、T2DM和silicosis是否形成具有公共卫生意义的重叠热点，并且糖尿病与矽肺这两类WHO关注的TB高风险状态是否在宿主髓系细胞中汇聚至可重复、可计算扰动的共同免疫调控程序？**
