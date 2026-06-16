# COVID-19 BALF scRNA-seq 分析流水线（Windows 适配版）

> 基于 Liao M. et al., *Single-cell landscape of bronchoalveolar immune cells in patients with COVID-19*, Nature Medicine, 2020.  
> 原始仓库: [zhangzlab/covid_balf](https://github.com/zhangzlab/covid_balf) | 在线数据: [covid19-balf.cells.ucsc.edu](https://covid19-balf.cells.ucsc.edu)

---

## 项目概述

本项目复现了上述论文的完整分析流程——从 66,452 个 BALF（支气管肺泡灌洗液）单细胞的聚类、注释，到巨噬细胞/NKT/T 细胞的亚群细分，再到细胞因子/趋化因子的跨疾病严重度比较。

**核心价值不在于"跑通了一个分析"，而在于从 Linux 服务器 → Windows 32GB 环境的整套适配与问题解决过程。**

---

## 思维过程：每个关键决策及其推理

### 决策 1：为什么不复用 Seurat 整合流程，而是直接用预计算的 integrated assay？

**原始脚本**（`seurat_integration.R`）在 Linux 服务器上对每个样本跑 `FindIntegrationAnchors` + `IntegrateData`，在巨噬细胞、NKT、T 细胞等子集上重复这一流程。

**在 Windows 32GB 下遇到的问题**：

```
Error in cholmod_allocate_factor: fatal error: out of memory
```

**问题诊断过程**：

1. 首先尝试 CCA 方法（`FindIntegrationAnchors(reduction = "cca")`）——CCA 在内部执行 `cbind` 合并稀疏矩阵时，Cholmod 库会请求一块**连续虚拟内存**。32GB Windows 下，物理内存 + 页文件有充足总空间，但没有单块足够大的连续地址空间。

2. 然后尝试 RPCA 方法（`reduction = "rpca"`）——需要对每个样本先跑 `ScaleData` + `RunPCA`。在大样本上成功，但在小样本（n < 200）导致 `AnnoySearch` 维度不匹配崩溃。为小样本排除添加了 `min_cells = 200` 过滤器。

3. 即使 RPCA 通过了 `FindIntegrationAnchors`，`IntegrateData` 的最后一步 `ScaleData` 仍然在合并后的大矩阵上触发 `std::bad_alloc`。这是物理内存不足，不是虚拟内存碎片问题。

**最终方案**：nCoV.rds 中已经有一个 `integrated` assay（作者预计算好的），包含 2000 个高变基因的整合矩阵。对于巨噬细胞/NKT/T 细胞的子集分析，直接复用这个 assay 做 PCA → 聚类，跳过重新整合：

```r
# 替代耗时的 integrate 流程：
DefaultAssay(Macrophage) <- "integrated"
Macrophage <- ScaleData(Macrophage, verbose = FALSE)
Macrophage <- RunPCA(Macrophage, verbose = FALSE)
Macrophage <- FindNeighbors(Macrophage, dims = 1:30)
Macrophage <- FindClusters(Macrophage, resolution = 0.8)
```

**代价评估**：这样做的降维空间更"全局"（2000 gene 是全局 HVG，不是巨噬细胞特异 HVG），但考虑到子集样本间的批次效应在全局整合中已被校正，这一折中在生物学上讲得通。对于已有整合数据的情况，这是合理的选择。

### 决策 2：MAST vs Wilcoxon 差异分析

**原始脚本**使用 MAST（`FindAllMarkers(test.use = "MAST")`）做差异分析。

**问题**：在 16 个巨噬细胞 cluster 上跑到第 7 个时发生 C++ segfault：

```
*** caught segfault *** address 0x7ff457bbe, cause 'memory not mapped'
```

MAST 是一个基于**跨细胞建模**的方法（hurdle model），每个 cluster 需要拟合一个 GLM，特别消耗内存在基因数 × 细胞数的全矩阵上。

**改成 `test.use = "wilcox"`（默认）**：Wilcoxon 秩和检验是非参数的、基于排名的，天然比 GLM 省内存。对于 >10K 细胞的数据，Wilcoxon 的统计效力与 MAST 高度相关，而且对 dropout 事件的容忍度更好。

**另一个坑**：Seurat v4 把列名从 `avg_logFC` 改成了 `avg_log2FC`。原脚本在引用差异结果时使用了旧列名，导致 `"找不到对象 'avg_logFC'"` 错误。

### 决策 3：动态 rename_map 模式

在 `run_NKT_fig.R` 和 `run_cytokine.R` 中，需要将数值 cluster ID 映射到生物学组名。直接硬编码一个 0-17 的映射表会遇到：

```
Cannot find identity 17
```

**原因**：如果某个 cluster 不存在于当前数据中（下采样、子集提取后 cluster 数量动态变化），`RenameIdents` 会崩溃。

**解决方案**——动态只映射实际存在的 cluster：

```r
mc <- levels(Idents(Macrophage.Integrated))
rename_map <- list()
for (cl in mc) {
  cl_num <- as.numeric(cl)
  if (cl_num == 1) rename_map[[cl]] <- 'Group1'
  else if (cl_num %in% c(13, 4, 8, 9, 0)) rename_map[[cl]] <- 'Group2'
  else if (cl_num %in% c(15, 5, 11)) rename_map[[cl]] <- 'Group3'
  else if (cl_num %in% c(2, 3, 6, 7, 10, 14)) rename_map[[cl]] <- 'Group4'
  else rename_map[[cl]] <- 'Other'
}
Macrophage.Integrated <- RenameIdents(Macrophage.Integrated, rename_map)
```

这个模式保证无论上游聚类产生了 12 个还是 18 个 cluster，映射都是完备的。

### 决策 4：顺序加载-计算-释放模式

`run_cytokine.R` 的第 2 部分需要比较 14 个趋化因子受体在巨噬细胞（4 组）、NKT（6 种）、其他 4 种免疫细胞上的平均表达。

**朴素做法**：同时加载 nCoV.rds（3.5GB） + 巨噬细胞 RDS + NKT RDS → 内存爆炸。

**采用的做法**：

```
加载 NKT RDS → AverageExpression() → 保存矩阵 → 释放
加载 nCoV.rds → 子集 Neutrophil/B/Plasma/pDC → AverageExpression() → 保存 → 释放
加载 巨噬细胞 RDS → AverageExpression() → 保存 → 释放
合并三个矩阵
```

每个 AverageExpression 只需要原始数据的对数归一化结果，不需要全部在内存中同时存在。这是典型的**时间换空间**策略。

---

## 流水线架构

```
nCoV.rds (66,452 cells × 22,000 genes)

┌─ scripts/adapted/run_total_fig.R ──── 全局细胞图谱
│   ├─ UMAP 可视化（全部细胞 / 按分组 / 按样本）
│   ├─ 12 个 marker 基因验证（FeaturePlot + Violin + DotPlot）
│   └─ 8 种免疫细胞比例统计（t-test: HC vs O vs S/C）
│
├─ scripts/adapted/run_macrophage.R ── 巨噬细胞亚群
│   ├─ 从 16 个 cluster 子集 → 下采样 → integrated assay 重聚类
│   └─ 输出: Myeloid/2-Macrophage.rds + marker 表格
│
├─ scripts/adapted/run_NKT_integration.R ── NK/T 细胞聚类
│   ├─ cluster 6,9,14,17 子集 → 下采样 → 14 个 subcluster
│   └─ 输出: NKT/6-NKT.rds
│
├─ scripts/adapted/run_NKT_fig.R ── NKT 注释与综合分析
│   ├─ 注释为 CD4 T / CD8 T / Cycling T / NK / Treg / innate T
│   ├─ 28 个 T/NK marker 的 DotPlot
│   ├─ 巨噬细胞 4 组分组: FCN1+/SPP1+/FABP4+/增殖型
│   └─ 输出: NKT/6-NKT-annotated.rds
│
├─ scripts/adapted/run_T_integration.R ── T 细胞精细亚群
│   ├─ 从 NKT 提取纯 T → 9 个 T 亚群
│   └─ 输出: NKT_T/6-T.rds
│
└─ scripts/adapted/run_cytokine.R ── 细胞因子/趋化因子分析
    ├─ Part 1: 22 种细胞因子在巨噬细胞 4 组上的表达
    ├─ Part 2: 14 种趋化因子受体跨 10 种细胞类型的比较
    └─ 三组间 DEG（Wilcoxon: S/C vs O vs HC）
```

## 文件结构

```
covid-balf-workflow/
├── README.md                          ← 你正在读的文件
├── CLAUDE.md                          ← Agent 工程技能配置
├── .gitignore
│
├── scripts/
│   ├── adapted/                       ← ★ Windows 适配版（本项目核心产出）
│   │   ├── run_total_fig.R            # 全局细胞图谱
│   │   ├── run_macrophage.R           # 巨噬细胞亚群
│   │   ├── run_NKT_integration.R      # NKT 聚类
│   │   ├── run_NKT_fig.R              # NKT 注释与可视化（最复杂）
│   │   ├── run_T_integration.R        # T 细胞亚群细分
│   │   ├── run_cytokine.R             # 细胞因子/趋化因子分析
│   │   ├── check_assays.R             # 诊断：检查 nCoV.rds 结构
│   │   ├── check_samples.R            # 诊断：检查样本分布
│   │   └── check_nkt.R                # 诊断：检查 NKT 子集
│   │
│   └── original/                      ← 原始 Zhang lab 脚本（保留出处）
│       ├── seurat_integration.R       # 全样本整合
│       ├── total_fig.R
│       ├── macrophage_fig.R
│       ├── macrophage_integration.R
│       ├── NKT.R
│       ├── NKT_integration.R
│       ├── T_integration.R
│       ├── cytokine_fig.R
│       └── tcr_fig.R
│
├── data/
│   ├── meta.txt                       # 样本 → 分组映射
│   ├── all.cell.annotation.meta.txt   # 全细胞类型注释
│   ├── myeloid.cell.annotation.meta.txt
│   ├── NKT.cell.annotation.meta.txt
│   └── README.original.md            # 原始 README
│
└── docs/agents/                       # 工程技能约定
    ├── issue-tracker.md
    ├── triage-labels.md
    └── domain.md
```

## 复现方法

### 1. 下载数据

```r
# 3.5GB，需要较长时间
download.file("http://cells.ucsc.edu/covid19-balf/nCoV.rds", "nCoV.rds")
```

### 2. 安装依赖

```r
install.packages(c("Seurat", "Matrix", "dplyr", "ggplot2", "reshape2", 
                   "ggpubr", "cowplot", "patchwork"))
```

### 3. 运行顺序

```r
source("scripts/adapted/run_total_fig.R")       # 1. 全局图谱
source("scripts/adapted/run_macrophage.R")      # 2. 巨噬细胞亚群
source("scripts/adapted/run_NKT_integration.R")  # 3. NKT 聚类
source("scripts/adapted/run_NKT_fig.R")         # 4. NKT 注释
source("scripts/adapted/run_T_integration.R")   # 5. T 细胞细分
source("scripts/adapted/run_cytokine.R")        # 6. 细胞因子分析
```

### 系统要求

- Windows 10/11 x64
- R 4.1.2+
- 矩阵物理内存 32GB+（建议）
- 如果内存不足，可减少 `run_macrophage.R` 中的 `max_cells` 参数

---

## 致谢

原始分析管道由 Zheng Zhang 实验室（深圳市第三人民医院/南方科技大学）开发。  
论文: Liao M. et al. *Nature Medicine* 26:842–844 (2020). [doi:10.1038/s41591-020-0901-9](https://doi.org/10.1038/s41591-020-0901-9)
