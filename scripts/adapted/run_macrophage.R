# Macrophage analysis adapted for Windows
# Based on macrophage_integration.R from Zhang lab
# Uses existing integrated assay from nCoV.rds to avoid memory overflow

.libPaths(c("E:/R/R_libs", .libPaths()))
library(Seurat)
library(Matrix)
library(dplyr)
library(ggplot2)
library(reshape2)

memory.limit(size = 32000)
gc()

# Paths
BASE_DIR <- "E:/Claude code/shengxin/covid_balf"
MYELOID_DIR <- file.path(BASE_DIR, "Myeloid")
dir.create(MYELOID_DIR, showWarnings = FALSE, recursive = TRUE)
setwd(MYELOID_DIR)

# Load pre-computed integrated object
nCoV.integrated <- readRDS(file.path(BASE_DIR, "nCoV.rds"))
cat("Loaded nCoV.rds:", dim(nCoV.integrated)[1], "genes x", dim(nCoV.integrated)[2], "cells\n")

# Subset macrophage clusters (clusters identified in the paper)
macrophage_clusters <- c('0','1','2','3','4','5','7','8','10','11','12','18','21','22','23','26')
Macrophage <- subset(nCoV.integrated, idents = macrophage_clusters)
cat("Macrophage subset:", dim(Macrophage)[1], "genes x", dim(Macrophage)[2], "cells\n")

write.table(table(Macrophage@meta.data$sample_new), file='Macrophage_sample.txt',
            quote = FALSE, sep='\t', row.names = FALSE)
write.table(table(Macrophage@meta.data$group), file='Macrophage_group.txt',
            quote = FALSE, sep='\t', row.names = FALSE)

# Downsample large samples for memory efficiency
max_cells_per_sample <- 2000
min_cells_per_sample <- 200
set.seed(42)
cells_to_keep <- c()
for (s in unique(Macrophage@meta.data$sample_new)) {
  cells_in_sample <- colnames(Macrophage)[Macrophage@meta.data$sample_new == s]
  nc <- length(cells_in_sample)
  if (nc < min_cells_per_sample) {
    cat(sprintf("Excluding %s: only %d cells (< %d)\n", s, nc, min_cells_per_sample))
    next
  }
  if (nc > max_cells_per_sample) {
    cells_to_keep <- c(cells_to_keep, sample(cells_in_sample, max_cells_per_sample))
    cat(sprintf("Downsampled %s: %d -> %d cells\n", s, nc, max_cells_per_sample))
  } else {
    cells_to_keep <- c(cells_to_keep, cells_in_sample)
  }
}
Macrophage <- subset(Macrophage, cells = cells_to_keep)
cat("After downsampling:", dim(Macrophage)[1], "genes x", dim(Macrophage)[2], "cells\n")
gc()

# Use existing integrated assay (from nCoV.rds) for clustering
# This avoids re-running expensive integration
DefaultAssay(Macrophage) <- "integrated"
Macrophage <- ScaleData(Macrophage, verbose = FALSE)
Macrophage <- RunPCA(Macrophage, verbose = FALSE)
Macrophage <- FindNeighbors(Macrophage, dims = 1:30)
Macrophage <- FindClusters(Macrophage, resolution = 0.8)
Macrophage <- RunTSNE(Macrophage, dims = 1:30)
Macrophage <- RunUMAP(Macrophage, reduction = "pca", dims = 1:30)
cat("Clustering done. Clusters:", length(unique(Idents(Macrophage))), "\n")
gc()

# Process RNA assay for marker finding
DefaultAssay(Macrophage) <- "RNA"
Macrophage[['percent.mito']] <- PercentageFeatureSet(Macrophage, pattern = "^MT-")
Macrophage <- NormalizeData(Macrophage, normalization.method = "LogNormalize", scale.factor = 1e4)
Macrophage <- FindVariableFeatures(Macrophage, selection.method = "vst", nfeatures = 2000, verbose = FALSE)
Macrophage <- ScaleData(Macrophage, verbose = FALSE)
cat("RNA processing done\n")
gc()

# QC plots
dpi <- 300
png(file="2-qc.png", width = dpi*16, height = dpi*8, units = "px", res = dpi)
print(patchwork::wrap_plots(VlnPlot(Macrophage, features = c("nFeature_RNA", "nCount_RNA"), combine = FALSE), ncol = 2))
dev.off()

png(file="2-umi-gene.png", width = dpi*6, height = dpi*5, units = "px", res = dpi)
print(FeatureScatter(Macrophage, feature1 = "nCount_RNA", feature2 = "nFeature_RNA"))
dev.off()

png(file="2-Macrophage-tsne.png", width = dpi*7, height = dpi*5, units = "px", res = dpi)
print(DimPlot(Macrophage, reduction = 'tsne', label = TRUE))
dev.off()

png(file="2-Macrophage-umap.png", width = dpi*7, height = dpi*5, units = "px", res = dpi)
print(DimPlot(Macrophage, reduction = 'umap', label = TRUE))
dev.off()

png(file="2-Macrophage-umap-split-sample.png", width = dpi*16, height = dpi*14, units = "px", res = dpi)
print(DimPlot(Macrophage, reduction = 'umap', label = TRUE, split.by = 'sample', ncol = 4))
dev.off()

png(file="2-Macrophage-umap-split-group.png", width = dpi*16, height = dpi*6, units = "px", res = dpi)
print(DimPlot(Macrophage, reduction = 'umap', label = TRUE, split.by = 'group', ncol = 3))
dev.off()

png(file="2-Macrophage-umap-group-sample.png", width = dpi*7, height = dpi*5, units = "px", res = dpi)
print(DimPlot(Macrophage, reduction = 'umap', label = TRUE, group.by = 'sample'))
dev.off()

png(file="2-Macrophage-umap-group-group.png", width = dpi*7, height = dpi*5, units = "px", res = dpi)
print(DimPlot(Macrophage, reduction = 'umap', label = TRUE, group.by = 'group'))
dev.off()

cat("All QC/UMAP plots saved\n")

# Save intermediate RDS before marker finding
saveRDS(Macrophage, file = "2-Macrophage.rds")
cat("Saved 2-Macrophage.rds\n")

# Find markers with Wilcoxon
DefaultAssay(Macrophage) <- "RNA"
cat("Finding macrophage markers...\n")
Macrophage@misc$markers <- FindAllMarkers(Macrophage, assay = 'RNA', only.pos = TRUE)
write.table(Macrophage@misc$markers, file='2-Macrophage-markers.txt',
            row.names = FALSE, quote = FALSE, sep='\t')
cat("Markers saved\n")

Macrophage@misc$averageExpression <- AverageExpression(Macrophage)
write.table(Macrophage@misc$averageExpression$RNA, file='2-Macrophage-average.txt',
            row.names = TRUE, quote = FALSE, sep='\t')

png(file="2-Macrophage-feature.png", width = dpi*24, height = dpi*5, units = "px", res = dpi)
print(VlnPlot(Macrophage, features = c("nFeature_RNA", "nCount_RNA")))
dev.off()

# Heatmap of top markers
hc.markers <- read.delim2("2-Macrophage-markers.txt", header = TRUE,
                          stringsAsFactors = FALSE, check.names = FALSE, sep = "\t")
hc.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC) -> top10
tt1 <- DoHeatmap(subset(Macrophage, downsample = 500), features = top10$gene) + NoLegend()
ggsave(file="2-Macrophage-feature2-1.pdf", plot = tt1, device = 'pdf',
       width = 20, height = 20, units = "in", dpi = dpi, limitsize = FALSE)

saveRDS(Macrophage, file = "2-Macrophage.rds")
cat("=== Macrophage analysis complete ===\n")
cat("Saved to:", MYELOID_DIR, "\n")
