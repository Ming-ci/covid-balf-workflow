# Overall cell landscape figures adapted for Windows
# Based on total_fig.R from Zhang lab (Liao et al., 2020, Nature Medicine)

.libPaths(c("E:/R/R_libs", .libPaths()))
library(Seurat)
library(Matrix)
library(dplyr)
library(ggplot2)
library(reshape2)
library(ggpubr)
library(cowplot)

BASE_DIR <- "E:/Claude code/shengxin/covid_balf"
FIGURE_DIR <- file.path(BASE_DIR, "Figure")
MARKER_DIR <- file.path(BASE_DIR, "marker")
dir.create(FIGURE_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(MARKER_DIR, showWarnings = FALSE, recursive = TRUE)
setwd(BASE_DIR)

nCoV.integrated <- readRDS(file.path(BASE_DIR, "nCoV.rds"))
cat("Loaded nCoV.rds:", dim(nCoV.integrated)[1], "genes x", dim(nCoV.integrated)[2], "cells\n")

dpi <- 300

# ==================== Fig 1a: UMAP of all cells ====================
cat("=== Fig 1a: UMAP of all cell clusters ===\n")
png(file = file.path(FIGURE_DIR, "1-umap_1.png"), width = dpi*8, height = dpi*6, units = "px", res = dpi)
pp <- DimPlot(nCoV.integrated, reduction = 'umap', pt.size = 0.5, label = FALSE, label.size = 5)
pp <- pp + theme(axis.title = element_text(size = 15), axis.text = element_text(size = 15, family = 'sans'),
                 legend.text = element_text(size = 15), axis.line = element_line(size = 0.8))
print(pp)
dev.off()
cat("  -> 1-umap_1.png\n")

# UMAP by group
png(file = file.path(FIGURE_DIR, "1-umap-by-group.png"), width = dpi*10, height = dpi*6, units = "px", res = dpi)
print(DimPlot(nCoV.integrated, reduction = 'umap', group.by = 'group', pt.size = 0.5))
dev.off()

# ==================== Cell type annotation ====================
cat("=== Cell type annotation ===\n")
nCoV.integrated <- RenameIdents(nCoV.integrated,
  '13'='Epithelial','16'='Epithelial','25'='Epithelial','28'='Epithelial','31'='Epithelial',
  '0'='Macrophages','1'='Macrophages','2'='Macrophages','3'='Macrophages','4'='Macrophages',
  '5'='Macrophages','7'='Macrophages','8'='Macrophages','10'='Macrophages','11'='Macrophages',
  '12'='Macrophages','18'='Macrophages','21'='Macrophages','22'='Macrophages','23'='Macrophages','26'='Macrophages',
  '30'='Mast','6'='T','9'='T','14'='T','17'='NK','19'='Plasma','27'='B',
  '15'='Neutrophil','20'='mDC','29'='pDC','24'='Doublets')
nCoV.integrated$celltype <- Idents(nCoV.integrated)
nCoV.integrated <- subset(nCoV.integrated, subset = celltype != 'Doublets')
nCoV.integrated$celltype <- factor(nCoV.integrated$celltype,
  levels = c('Epithelial','Macrophages','Neutrophil','mDC','pDC','Mast','T','NK','B','Plasma'))
Idents(nCoV.integrated) <- nCoV.integrated$celltype

# ==================== UMAP split by group ====================
cat("=== UMAP split by group ===\n")
png(file = file.path(FIGURE_DIR, "1-umap-split-group.png"), width = dpi*11, height = dpi*4.3, units = "px", res = dpi)
print(DimPlot(nCoV.integrated, reduction = 'umap', label = FALSE, split.by = 'group', ncol = 3, repel = TRUE))
dev.off()

# ==================== UMAP split by sample ====================
cat("=== UMAP split by sample ===\n")
png(file = file.path(FIGURE_DIR, "1-umap-split-sample.png"), width = dpi*18, height = dpi*9, units = "px", res = dpi)
print(DimPlot(nCoV.integrated, reduction = 'umap', label = FALSE, split.by = 'sample_new', ncol = 5))
dev.off()

# ==================== Marker gene feature plots ====================
cat("=== Marker gene plots ===\n")
markers <- c('TPPP3','KRT18','CD68','FCGR3B','CD1C','CLEC9A','LILRA4','TPSB2','CD3D','KLRD1','MS4A1','IGHG4')

# UMAP feature plots
existing_markers <- markers[markers %in% rownames(nCoV.integrated)]
cat("Plotting", length(existing_markers), "markers\n")
fp_list <- FeaturePlot(nCoV.integrated, features = existing_markers, cols = c("lightgrey","#ff0000"), combine = FALSE)
png(file = file.path(MARKER_DIR, "umap_markers.png"), width = dpi*16, height = dpi*10, units = "px", res = dpi)
print(patchwork::wrap_plots(fp_list, ncol = 4))
dev.off()

# Violin plots
vp_list <- VlnPlot(nCoV.integrated, features = existing_markers, pt.size = 0, combine = FALSE)
vp_list <- lapply(vp_list, function(p) p + theme(axis.text.x = element_text(angle = 45, hjust = 1)))
png(file = file.path(MARKER_DIR, "violin_markers.png"), width = dpi*24, height = dpi*18, units = "px", res = dpi)
print(patchwork::wrap_plots(vp_list, ncol = 4))
dev.off()

# Dot plot
pdf(file = file.path(MARKER_DIR, "dotplot_markers.pdf"), width = 10, height = 6)
pp <- DotPlot(nCoV.integrated, features = rev(markers), cols = c('white','#F8766D'), dot.scale = 5) + RotatedAxis()
pp <- pp + theme(axis.text.x = element_text(size = 12), axis.text.y = element_text(size = 12)) +
  labs(x = '', y = '') +
  guides(color = guide_colorbar(title = 'Scale expression'), size = guide_legend(title = 'Percent expressed'))
print(pp)
dev.off()

# ==================== Cell proportion analysis ====================
cat("=== Cell proportion analysis ===\n")
# Create sample-to-group mapping by merging sample_new with group info
meta <- nCoV.integrated@meta.data
sample_group_map <- unique(meta[, c("sample_new", "group")])
rownames(sample_group_map) <- NULL

# Cell type percentage per sample
cell_per_sample <- as.data.frame.matrix(table(meta$sample_new, meta$celltype))
cell_per_sample <- cell_per_sample / rowSums(cell_per_sample)
cell_per_sample$sample <- rownames(cell_per_sample)
cell_per_sample <- left_join(cell_per_sample, sample_group_map, by = c("sample" = "sample_new"))

# Proportion bar plot per sample
cell_per_sample_long <- melt(cell_per_sample, id.vars = c("sample", "group"),
                              measure.vars = setdiff(colnames(cell_per_sample), c("sample", "group")))
colnames(cell_per_sample_long) <- c("sample", "group", "celltype", "proportion")
cell_per_sample_long$proportion <- as.numeric(cell_per_sample_long$proportion)

# Bar plot of cell composition by sample
samples_ordered <- c('HC1','HC2','HC3','HC4','O1','O2','O3','S1','C1','C2','C3','C4','C5')
cell_per_sample_long$sample <- factor(cell_per_sample_long$sample, levels = samples_ordered)

cols_celltype <- c('#32b8ec','#60c3f0','#8ccdf1','#cae5f7','#92519c','#b878b0','#d7b1d2',
                   '#e7262a','#e94746','#eb666d','#ee838f','#f4abac','#fad9d9')

png(file = file.path(FIGURE_DIR, "1-cell-proportion-by-sample.png"), width = dpi*10, height = dpi*5, units = "px", res = dpi)
pp <- ggplot(cell_per_sample_long, aes(x = sample, y = proportion, fill = celltype)) +
  geom_bar(stat = "identity", width = 0.6, position = position_fill(reverse = TRUE), size = 0.2, colour = '#222222') +
  labs(x = '', y = 'Fraction of cells') +
  theme_cowplot() +
  theme(axis.text.x = element_text(size = 10, angle = 45, hjust = 1),
        axis.text.y = element_text(size = 12),
        legend.text = element_text(size = 10),
        legend.key.height = unit(4, 'mm')) +
  scale_fill_brewer(palette = "Set3")
print(pp)
dev.off()

# ==================== Statistical comparison of cell proportions ====================
cat("=== Statistical comparison ===\n")
cell_per_group <- as.data.frame.matrix(table(meta$sample_new, meta$celltype))
cell_per_group <- cell_per_group / rowSums(cell_per_group)
cell_per_group$sample <- rownames(cell_per_group)
cell_per_group <- left_join(cell_per_group, sample_group_map, by = c("sample" = "sample_new"))

immune_types <- c('Macrophages','Neutrophil','mDC','pDC','T','NK','B','Plasma')
pplist <- list()

for (ct in immune_types) {
  ct_data <- cell_per_group[, c("sample", "group", ct)]
  colnames(ct_data)[3] <- "percent"
  ct_data$percent <- as.numeric(ct_data$percent)

  ct_data <- ct_data %>% group_by(group) %>%
    mutate(upper = quantile(percent, 0.75, na.rm = TRUE),
           lower = quantile(percent, 0.25, na.rm = TRUE),
           mean = mean(percent, na.rm = TRUE),
           median = median(percent, na.rm = TRUE))

  my_comparisons <- list(c("HC", "O"), c("O", "S/C"), c("HC", "S/C"))

  pp1 <- ggplot(ct_data, aes(x = group, y = percent)) +
    geom_jitter(shape = 21, aes(fill = group), width = 0.25) +
    stat_summary(fun = mean, geom = "point", color = "grey60") +
    theme_cowplot() +
    theme(axis.text = element_text(size = 8), axis.title = element_text(size = 8),
          legend.position = 'none') +
    labs(title = ct, y = 'Percentage') +
    geom_errorbar(aes(ymin = lower, ymax = upper), col = "grey60", width = 0.25) +
    stat_compare_means(comparisons = my_comparisons, size = 2.5, method = "t.test")

  pplist[[ct]] <- pp1
}

pdf(file = file.path(FIGURE_DIR, "cell-percentage-stats.pdf"), width = 10, height = 5)
print(plot_grid(pplist[['Macrophages']], pplist[['Neutrophil']],
                pplist[['mDC']], pplist[['pDC']],
                pplist[['T']], pplist[['NK']], pplist[['B']],
                pplist[['Plasma']], ncol = 4, nrow = 2))
dev.off()

cat("\n=== Figure generation complete ===\n")
cat("Figures saved to:", FIGURE_DIR, "\n")
cat("Markers saved to:", MARKER_DIR, "\n")
