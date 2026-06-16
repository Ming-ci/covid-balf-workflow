# NKT cell type annotation and figure generation adapted for Windows
# Based on NKT.R (NKT_fig) from Zhang lab

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
FIGURE_DIR <- file.path(BASE_DIR, "Figure")
NKT_DIR <- file.path(BASE_DIR, "NKT")
dir.create(FIGURE_DIR, showWarnings = FALSE, recursive = TRUE)
setwd(FIGURE_DIR)

dpi <- 300

# ============ Load NKT data ============
NKT.Integrated <- readRDS(file.path(NKT_DIR, "6-NKT.rds"))
cat("Loaded 6-NKT.rds:", dim(NKT.Integrated)[1], "genes x", dim(NKT.Integrated)[2], "cells\n")
cat("Clusters:", paste(levels(Idents(NKT.Integrated)), collapse=", "), "\n")

# ============ Cell type annotation ============
# Rename clusters based on original paper's annotation
# Original mapping: 0=CCR7+ T, 1=CD8 T, 2=CD8 T, 3=Cycling T, 4=NK, 5=CD8 T,
#                   6=Treg, 7=Doublets, 8=NK, 9=innate T, 10=Cycling T, 11-12=Uncertain, 13=Doublets
# Our clustering may differ slightly - adapt based on our 14 clusters
NKT.Integrated <- RenameIdents(object = NKT.Integrated,
  '0' = 'CD4 T', '1' = 'CD8 T', '2' = 'CD8 T', '3' = 'Cycling T',
  '4' = 'NK', '5' = 'CD8 T', '6' = 'Treg', '7' = 'NK',
  '8' = 'innate T', '9' = 'Cycling T',
  '10' = 'Uncertain', '11' = 'Uncertain', '12' = 'Uncertain', '13' = 'Uncertain')
NKT.Integrated$celltype <- Idents(NKT.Integrated)
nkt_celltypes <- c('CD4 T', 'CD8 T', 'Cycling T', 'NK', 'Treg', 'innate T', 'Uncertain')
NKT.Integrated$celltype <- factor(NKT.Integrated$celltype, levels = nkt_celltypes)
Idents(NKT.Integrated) <- NKT.Integrated$celltype

# Remove Uncertain
NKT.clean <- subset(NKT.Integrated, idents = c('CD4 T', 'CD8 T', 'Cycling T', 'NK', 'Treg', 'innate T'))
cat("NKT clean (no Uncertain):", ncol(NKT.clean), "cells in", length(unique(Idents(NKT.clean))), "types\n")

# ============ UMAP of NKT cell types ============
cols_nkt <- c('#f8766d', '#cd9600', '#7cae00', '#00be67', '#00bfc4', '#00a9ff')

png(file = "4-umap_1.png", width = dpi*6, height = dpi*4, units = "px", res = dpi)
pp <- DimPlot(object = NKT.clean, reduction = 'umap', label = FALSE, pt.size = 0.8,
              cols = cols_nkt, repel = TRUE)
pp <- pp + theme(axis.title = element_text(size = 12), axis.text = element_text(size = 12),
                  legend.text = element_text(size = 12))
print(pp)
dev.off()
cat("  -> 4-umap_1.png\n")

# ============ UMAP split by group ============
png(file = "4-NKT-umap-split-group.png", width = dpi*11, height = dpi*4.3, units = "px", res = dpi)
pp_temp <- DimPlot(object = NKT.clean, reduction = 'umap', label = FALSE,
                    split.by = 'group', ncol = 3, repel = TRUE, pt.size = 1.2,
                    cols = cols_nkt)
pp_temp <- pp_temp + theme(axis.title = element_text(size = 17), axis.text = element_text(size = 17),
                            legend.text = element_text(size = 17))
print(pp_temp)
dev.off()
cat("  -> 4-NKT-umap-split-group.png\n")

# ============ Marker gene dotplot ============
nCoV_markers <- c('CD3D','IL7R','CCR7','GZMA','CD8A','CD8B','CXCR3','GZMK',
                   'MKI67','TYMS','NKG7','GZMB','KLRD1','KLRC1','XCL1','KLRF1',
                   'EOMES','CX3CR1','CD4','FOXP3','CTLA4','IL2RA','CXCR6',
                   'TRGV9','TRDV2','KLRB1','SLC4A10')
existing_markers <- nCoV_markers[nCoV_markers %in% rownames(NKT.clean)]
cat("Plotting", length(existing_markers), "/", length(nCoV_markers), "NKT markers\n")

pdf(file = "4-heatmap-seurat.pdf", width = 8.5, height = 3.5)
pp <- DotPlot(NKT.clean, features = rev(existing_markers), cols = c('white', '#F8766D'), dot.scale = 5) + RotatedAxis()
pp <- pp + theme(axis.text.x = element_text(size = 10), axis.text.y = element_text(size = 12)) +
  labs(x = '', y = '') +
  guides(color = guide_colorbar(title = 'Scale expression'), size = guide_legend(title = 'Percent expressed'))
print(pp)
dev.off()
cat("  -> 4-heatmap-seurat.pdf\n")

# ============ Cell proportion per sample ============
NKT.clean[["cluster"]] <- Idents(NKT.clean)
cell_summary <- table(NKT.clean@meta.data$sample_new, NKT.clean@meta.data$cluster)
write.table(cell_summary, file = '4-NKT-percentage-sample.txt', quote = FALSE, sep = '\t')

# Bar plot of cell proportions by sample
organ.summary <- as.data.frame.matrix(cell_summary)
organ.summary$sample <- rownames(organ.summary)
organ.summary.long <- melt(organ.summary, id.vars = "sample")
colnames(organ.summary.long) <- c("sample", "celltype", "count")
organ.summary.long$count <- as.numeric(organ.summary.long$count)

# Add group info
meta <- NKT.clean@meta.data
sample_group_map <- unique(meta[, c("sample_new", "group")])
organ.summary.long <- left_join(organ.summary.long, sample_group_map,
                                 by = c("sample" = "sample_new"))

pdf(file = "4-NKT-sample-percentage.pdf", width = 7, height = 3.6)
pp <- ggplot(data = organ.summary.long, aes(x = sample, y = count, fill = celltype)) +
  geom_bar(stat = "identity", width = 0.6, position = position_fill(reverse = FALSE),
           size = 0.5, colour = '#222222') +
  labs(x = '', y = 'Fraction of cells') +
  cowplot::theme_cowplot() +
  theme(axis.text.x = element_text(angle = 45, size = 10, hjust = 0.5, vjust = 0.5),
        axis.text.y = element_text(size = 14),
        legend.text = element_text(size = 14),
        legend.key.height = unit(5, 'mm')) +
  scale_fill_manual(values = cols_nkt)
print(pp)
dev.off()
cat("  -> 4-NKT-sample-percentage.pdf\n")

# ============ Find markers per annotated cell type ============
DefaultAssay(NKT.clean) <- "RNA"
cat("Finding markers for NKT cell types...\n")
NKT.clean@misc$markers <- FindAllMarkers(NKT.clean, assay = 'RNA', only.pos = TRUE)
write.table(NKT.clean@misc$markers, file = '4-markers.txt',
            row.names = FALSE, quote = FALSE, sep = '\t')
cat("  -> 4-markers.txt\n")

# ============ Heatmap of top markers ============
hc.markers <- read.delim2("4-markers.txt", header = TRUE,
                          stringsAsFactors = FALSE, check.names = FALSE, sep = "\t")
nkt_groups_rev <- c('CD4 T', 'CD8 T', 'Cycling T', 'NK', 'Treg', 'innate T')
hc.markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_log2FC) -> top20
# Deduplicate genes
hvg <- unique(top20$gene)
cat("Top marker genes for heatmap:", length(hvg), "\n")

var.genes <- unique(c(NKT.clean@assays$RNA@var.features, top20$gene))
NKT.clean <- ScaleData(NKT.clean, verbose = FALSE, features = var.genes)
tt1 <- DoHeatmap(object = NKT.clean, features = hvg, angle = 0, hjust = 0.5, size = 6) + NoLegend()
ggsave(file = "4-feature2-1.pdf", plot = tt1, device = 'pdf',
       width = 20, height = 24, units = "in", dpi = dpi, limitsize = FALSE)
cat("  -> 4-feature2-1.pdf\n")

# Save annotated NKT before freeing memory
Idents(NKT.clean) <- NKT.clean$celltype
saveRDS(NKT.clean, file = file.path(NKT_DIR, "6-NKT-annotated.rds"))
cat("Saved 6-NKT-annotated.rds\n")

# ============ Macrophage vs NKT comparison ============
# Free memory before loading more data
rm(NKT.Integrated, NKT.clean)
gc()
cat("=== Macrophage vs NKT comparison ===\n")

# Annotate macrophage subtypes (based on paper grouping)
# Our macrophage data has 16 clusters (0-15), adapt annotation accordingly
Macrophage.Integrated <- readRDS(file.path(BASE_DIR, "Myeloid", "2-Macrophage.rds"))
cat("Loaded macrophage RDS:", ncol(Macrophage.Integrated), "cells, clusters:", paste(levels(Idents(Macrophage.Integrated)), collapse=", "), "\n")

# Check which clusters exist
mc <- levels(Idents(Macrophage.Integrated))
cat("Existing clusters:", paste(mc, collapse=", "), "\n")

# Build rename mapping from existing clusters
rename_map <- list()
for (cl in mc) {
  cl_num <- as.numeric(cl)
  if (cl_num == 1) rename_map[[cl]] <- 'Group1'
  else if (cl_num %in% c(13, 4, 8, 9, 0)) rename_map[[cl]] <- 'Group2'
  else if (cl_num %in% c(15, 5, 11)) rename_map[[cl]] <- 'Group3'
  else if (cl_num %in% c(2, 3, 6, 7, 10, 14)) rename_map[[cl]] <- 'Group4'
  else rename_map[[cl]] <- 'Other'
}
cat("Rename map:\n")
for (k in names(rename_map)) cat("  Cluster", k, "->", rename_map[[k]], "\n")

Macrophage.Integrated <- RenameIdents(Macrophage.Integrated, rename_map)
Macrophage.Integrated$macro_group <- Idents(Macrophage.Integrated)
Macrophage.Integrated <- subset(Macrophage.Integrated, subset = macro_group != 'Other')

# Reload NKT data
NKT.clean <- readRDS(file.path(NKT_DIR, "6-NKT-annotated.rds"))
# In case the annotated version doesn't exist, use the original
if (!"celltype" %in% colnames(NKT.clean@meta.data)) {
  NKT.clean <- RenameIdents(NKT.clean,
    '0' = 'CD4 T', '1' = 'CD8 T', '2' = 'CD8 T', '3' = 'Cycling T',
    '4' = 'NK', '5' = 'CD8 T', '6' = 'Treg', '7' = 'NK',
    '8' = 'innate T', '9' = 'Cycling T',
    '10' = 'Uncertain', '11' = 'Uncertain', '12' = 'Uncertain', '13' = 'Uncertain')
  NKT.clean$celltype <- Idents(NKT.clean)
  NKT.clean <- subset(NKT.clean, idents = c('CD4 T', 'CD8 T', 'Cycling T', 'NK', 'Treg', 'innate T'))
}

cat("NKT clean cells:", ncol(NKT.clean), "\n")

# Simple cell proportion comparison between groups
# Compare MP vs T proportions in Severe vs Moderate
macro_meta <- Macrophage.Integrated@meta.data
nkt_meta <- NKT.clean@meta.data

# For this comparison, use the group annotation from metadata
cat("Severe (S/C) vs Moderate (O) comparison:\n")
cat("  Macrophage cells in S/C:", sum(macro_meta$group == 'S/C'), "\n")
cat("  Macrophage cells in O:", sum(macro_meta$group == 'O'), "\n")
cat("  NKT cells in S/C:", sum(nkt_meta$group == 'S/C'), "\n")
cat("  NKT cells in O:", sum(nkt_meta$group == 'O'), "\n")

write.table(macro_meta[, c("sample_new", "group", "macro_group")],
            file = '4-macro-group-metadata.txt', row.names = FALSE, quote = FALSE, sep = '\t')
write.table(nkt_meta[, c("sample_new", "group", "celltype")],
            file = '4-nkt-group-metadata.txt', row.names = FALSE, quote = FALSE, sep = '\t')
cat("  -> 4-macro-group-metadata.txt, 4-nkt-group-metadata.txt\n")

# ============ Cytokine/cytof markers dotplot ============
cytokines <- c('CD33','CD14','CXCR6','CD163','IL3RA','CD27','CD19','CCR6',
               'NCAM1','KLRB1','CD69','CD68','CXCR3','PDCD1','CCR7','ITGAE',
               'CX3CR1','CD4','FCGR3A','IL7R','ITGAM','ITGAX','CD38','PTPRC',
               'B3GAT1','CD3D','CD8A')
existing_cyto <- cytokines[cytokines %in% rownames(NKT.clean)]

pdf(file = "4-heatmap-cytokine_cytof_sample.pdf", width = 10, height = 4)
pp <- DotPlot(NKT.clean, features = rev(existing_cyto), cols = c('white', '#F8766D'), dot.scale = 5) + RotatedAxis()
pp <- pp + theme(axis.text.x = element_text(size = 9), axis.text.y = element_text(size = 9)) +
  labs(x = '', y = '') +
  guides(color = guide_colorbar(title = 'Scale expression'), size = guide_legend(title = 'Percent expressed'))
print(pp)
dev.off()
cat("  -> 4-heatmap-cytokine_cytof_sample.pdf\n")

# ============ Final save ============
# Annotated NKT already saved above; just save the final updated version
saveRDS(NKT.clean, file = file.path(NKT_DIR, "6-NKT-annotated.rds"))

cat("\n=== NKT figure generation complete ===\n")
cat("Figures saved to:", FIGURE_DIR, "\n")
