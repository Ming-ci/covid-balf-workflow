# T cell subset analysis adapted for Windows
# Based on T_integration.R from Zhang lab
# Subsets T cells from NKT and re-clusters

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
NKT_T_DIR <- file.path(BASE_DIR, "NKT_T")
NKT_DIR <- file.path(BASE_DIR, "NKT")
dir.create(NKT_T_DIR, showWarnings = FALSE, recursive = TRUE)
setwd(NKT_T_DIR)

# Load annotated NKT data and subset T cells
NKT <- readRDS(file.path(NKT_DIR, "6-NKT-annotated.rds"))
cat("Loaded NKT:", dim(NKT)[1], "genes x", dim(NKT)[2], "cells\n")

# Check cell types
if ("celltype" %in% colnames(NKT@meta.data)) {
  Idents(NKT) <- NKT$celltype
  cat("Cell types:", paste(unique(NKT$celltype), collapse=", "), "\n")
  T_cell_types <- c('CD4 T', 'CD8 T', 'Cycling T', 'Treg', 'innate T')
  T_cells <- subset(NKT, idents = T_cell_types)
} else {
  # Fallback: subset known T cell clusters from original NKT
  T_cells <- subset(NKT, idents = c('0','1','2','3','5','6','8','9'))
}
cat("T cell subset:", dim(T_cells)[1], "genes x", dim(T_cells)[2], "cells\n")

# Free NKT from memory
rm(NKT)
gc()

write.table(table(T_cells@meta.data$sample_new), file='T_sample.txt',
            quote = FALSE, sep='\t', row.names = FALSE)
write.table(table(T_cells@meta.data$group), file='T_group.txt',
            quote = FALSE, sep='\t', row.names = FALSE)

# Downsample if needed
max_cells_per_sample <- 2000
min_cells_per_sample <- 30
set.seed(42)
cells_to_keep <- c()
for (s in unique(T_cells@meta.data$sample_new)) {
  cells_in_sample <- colnames(T_cells)[T_cells@meta.data$sample_new == s]
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
T_cells <- subset(T_cells, cells = cells_to_keep)
cat("After downsampling:", dim(T_cells)[1], "genes x", dim(T_cells)[2], "cells\n")
gc()

# Use existing integrated assay for clustering
DefaultAssay(T_cells) <- "integrated"
T_cells <- ScaleData(T_cells, verbose = FALSE)
T_cells <- RunPCA(T_cells, verbose = FALSE)
T_cells <- FindNeighbors(T_cells, dims = 1:30)
T_cells <- FindClusters(T_cells, resolution = 0.8)
T_cells <- RunTSNE(T_cells, dims = 1:30)
T_cells <- RunUMAP(T_cells, reduction = "pca", dims = 1:30)
cat("Clustering done. Clusters:", length(unique(Idents(T_cells))), "\n")
gc()

# Process RNA assay
DefaultAssay(T_cells) <- "RNA"
T_cells[['percent.mito']] <- PercentageFeatureSet(T_cells, pattern = "^MT-")
T_cells <- NormalizeData(T_cells, normalization.method = "LogNormalize", scale.factor = 1e4)
T_cells <- FindVariableFeatures(T_cells, selection.method = "vst", nfeatures = 2000, verbose = FALSE)
T_cells <- ScaleData(T_cells, verbose = FALSE)
cat("RNA processing done\n")
gc()

# QC plots
dpi <- 300
png(file="T-qc.png", width = dpi*16, height = dpi*8, units = "px", res = dpi)
print(patchwork::wrap_plots(VlnPlot(T_cells, features = c("nFeature_RNA", "nCount_RNA"), combine = FALSE), ncol = 2))
dev.off()

png(file="T-umi-gene.png", width = dpi*6, height = dpi*5, units = "px", res = dpi)
print(FeatureScatter(T_cells, feature1 = "nCount_RNA", feature2 = "nFeature_RNA"))
dev.off()

png(file="T-tsne.png", width = dpi*7, height = dpi*5, units = "px", res = dpi)
print(DimPlot(T_cells, reduction = 'tsne', label = TRUE))
dev.off()

png(file="T-umap.png", width = dpi*7, height = dpi*5, units = "px", res = dpi)
print(DimPlot(T_cells, reduction = 'umap', label = TRUE))
dev.off()

png(file="T-umap-split-sample.png", width = dpi*16, height = dpi*14, units = "px", res = dpi)
print(DimPlot(T_cells, reduction = 'umap', label = TRUE, split.by = 'sample', ncol = 4))
dev.off()

png(file="T-umap-split-group.png", width = dpi*16, height = dpi*6, units = "px", res = dpi)
print(DimPlot(T_cells, reduction = 'umap', label = TRUE, split.by = 'group', ncol = 3))
dev.off()

png(file="T-umap-group-sample.png", width = dpi*7, height = dpi*5, units = "px", res = dpi)
print(DimPlot(T_cells, reduction = 'umap', label = TRUE, group.by = 'sample'))
dev.off()

png(file="T-umap-group-group.png", width = dpi*7, height = dpi*5, units = "px", res = dpi)
print(DimPlot(T_cells, reduction = 'umap', label = TRUE, group.by = 'group'))
dev.off()

cat("All QC/UMAP plots saved\n")

# Save intermediate
saveRDS(T_cells, file = "6-T.rds")
cat("Saved 6-T.rds\n")

# Find markers
DefaultAssay(T_cells) <- "RNA"
cat("Finding T cell markers...\n")
T_cells@misc$markers <- FindAllMarkers(T_cells, assay = 'RNA', only.pos = TRUE)
write.table(T_cells@misc$markers, file='T-markers.txt',
            row.names = FALSE, quote = FALSE, sep='\t')
cat("Markers saved\n")

T_cells@misc$averageExpression <- AverageExpression(T_cells)
write.table(T_cells@misc$averageExpression$RNA, file='T-average.txt',
            row.names = TRUE, quote = FALSE, sep='\t')

png(file="T-feature.png", width = dpi*24, height = dpi*5, units = "px", res = dpi)
print(VlnPlot(T_cells, features = c("nFeature_RNA", "nCount_RNA")))
dev.off()

# Heatmap
hc.markers <- read.delim2("T-markers.txt", header = TRUE,
                          stringsAsFactors = FALSE, check.names = FALSE, sep = "\t")
hc.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC) -> top10
tt1 <- DoHeatmap(subset(T_cells, downsample = 500), features = top10$gene) + NoLegend()
ggsave(file="T-feature2-1.pdf", plot = tt1, device = 'pdf',
       width = 20, height = 16, units = "in", dpi = dpi, limitsize = FALSE)

saveRDS(T_cells, file = "6-T.rds")
cat("=== T cell analysis complete ===\n")
cat("Saved to:", NKT_T_DIR, "\n")
