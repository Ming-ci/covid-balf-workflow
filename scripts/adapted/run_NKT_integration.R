# NKT integration analysis adapted for Windows
# Based on NKT_integration.R from Zhang lab
# Uses existing integrated assay from nCoV.rds

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
NKT_DIR <- file.path(BASE_DIR, "NKT")
dir.create(NKT_DIR, showWarnings = FALSE, recursive = TRUE)
setwd(NKT_DIR)

# Load pre-computed integrated object
nCoV.integrated <- readRDS(file.path(BASE_DIR, "nCoV.rds"))
cat("Loaded nCoV.rds:", dim(nCoV.integrated)[1], "genes x", dim(nCoV.integrated)[2], "cells\n")

# Subset NKT clusters (clusters 6,9,14,17 from original analysis = T+NK)
NKT <- subset(nCoV.integrated, idents = c('6','9','14','17'))
cat("NKT subset:", dim(NKT)[1], "genes x", dim(NKT)[2], "cells\n")

write.table(table(NKT@meta.data$sample_new), file='NKT_sample.txt',
            quote = FALSE, sep='\t', row.names = FALSE)
write.table(table(NKT@meta.data$group), file='NKT_group.txt',
            quote = FALSE, sep='\t', row.names = FALSE)

# Downsample large samples for memory efficiency
max_cells_per_sample <- 2000
min_cells_per_sample <- 50
set.seed(42)
cells_to_keep <- c()
for (s in unique(NKT@meta.data$sample_new)) {
  cells_in_sample <- colnames(NKT)[NKT@meta.data$sample_new == s]
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
NKT <- subset(NKT, cells = cells_to_keep)
cat("After downsampling:", dim(NKT)[1], "genes x", dim(NKT)[2], "cells\n")
gc()

# Use existing integrated assay for clustering
DefaultAssay(NKT) <- "integrated"
NKT <- ScaleData(NKT, verbose = FALSE)
NKT <- RunPCA(NKT, verbose = FALSE)
NKT <- FindNeighbors(NKT, dims = 1:30)
NKT <- FindClusters(NKT, resolution = 0.8)
NKT <- RunTSNE(NKT, dims = 1:30)
NKT <- RunUMAP(NKT, reduction = "pca", dims = 1:30)
cat("Clustering done. Clusters:", length(unique(Idents(NKT))), "\n")
gc()

# Process RNA assay for marker finding
DefaultAssay(NKT) <- "RNA"
NKT[['percent.mito']] <- PercentageFeatureSet(NKT, pattern = "^MT-")
NKT <- NormalizeData(NKT, normalization.method = "LogNormalize", scale.factor = 1e4)
NKT <- FindVariableFeatures(NKT, selection.method = "vst", nfeatures = 2000, verbose = FALSE)
NKT <- ScaleData(NKT, verbose = FALSE)
cat("RNA processing done\n")
gc()

# QC plots
dpi <- 300
png(file="6-qc.png", width = dpi*16, height = dpi*8, units = "px", res = dpi)
print(patchwork::wrap_plots(VlnPlot(NKT, features = c("nFeature_RNA", "nCount_RNA"), combine = FALSE), ncol = 2))
dev.off()

png(file="6-umi-gene.png", width = dpi*6, height = dpi*5, units = "px", res = dpi)
print(FeatureScatter(NKT, feature1 = "nCount_RNA", feature2 = "nFeature_RNA"))
dev.off()

png(file="6-NKT-tsne.png", width = dpi*7, height = dpi*5, units = "px", res = dpi)
print(DimPlot(NKT, reduction = 'tsne', label = TRUE))
dev.off()

png(file="6-NKT-umap.png", width = dpi*7, height = dpi*5, units = "px", res = dpi)
print(DimPlot(NKT, reduction = 'umap', label = TRUE))
dev.off()

png(file="6-NKT-umap-split-sample.png", width = dpi*16, height = dpi*14, units = "px", res = dpi)
print(DimPlot(NKT, reduction = 'umap', label = TRUE, split.by = 'sample', ncol = 4))
dev.off()

png(file="6-NKT-umap-split-group.png", width = dpi*16, height = dpi*6, units = "px", res = dpi)
print(DimPlot(NKT, reduction = 'umap', label = TRUE, split.by = 'group', ncol = 3))
dev.off()

png(file="6-NKT-umap-group-sample.png", width = dpi*7, height = dpi*5, units = "px", res = dpi)
print(DimPlot(NKT, reduction = 'umap', label = TRUE, group.by = 'sample'))
dev.off()

png(file="6-NKT-umap-group-group.png", width = dpi*7, height = dpi*5, units = "px", res = dpi)
print(DimPlot(NKT, reduction = 'umap', label = TRUE, group.by = 'group'))
dev.off()

cat("All QC/UMAP plots saved\n")

# Save intermediate RDS before marker finding
saveRDS(NKT, file = "6-NKT.rds")
cat("Saved 6-NKT.rds\n")

# Find markers with Wilcoxon
DefaultAssay(NKT) <- "RNA"
cat("Finding NKT markers...\n")
NKT@misc$markers <- FindAllMarkers(NKT, assay = 'RNA', only.pos = TRUE)
write.table(NKT@misc$markers, file='6-NKT-markers.txt',
            row.names = FALSE, quote = FALSE, sep='\t')
cat("Markers saved\n")

NKT@misc$averageExpression <- AverageExpression(NKT)
write.table(NKT@misc$averageExpression$RNA, file='6-NKT-average.txt',
            row.names = TRUE, quote = FALSE, sep='\t')

png(file="6-NKT-feature.png", width = dpi*24, height = dpi*5, units = "px", res = dpi)
print(VlnPlot(NKT, features = c("nFeature_RNA", "nCount_RNA")))
dev.off()

# Heatmap of top markers
hc.markers <- read.delim2("6-NKT-markers.txt", header = TRUE,
                          stringsAsFactors = FALSE, check.names = FALSE, sep = "\t")
hc.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC) -> top10
tt1 <- DoHeatmap(subset(NKT, downsample = 500), features = top10$gene) + NoLegend()
ggsave(file="6-NKT-feature2-1.pdf", plot = tt1, device = 'pdf',
       width = 20, height = 10, units = "in", dpi = dpi, limitsize = FALSE)

saveRDS(NKT, file = "6-NKT.rds")
cat("=== NKT analysis complete ===\n")
cat("Saved to:", NKT_DIR, "\n")
