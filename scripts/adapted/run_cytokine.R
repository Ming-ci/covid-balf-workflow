# Cytokine analysis adapted for Windows
# Based on cytokine_fig.R from Zhang lab

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
MYELOID_DIR <- file.path(BASE_DIR, "Myeloid")
NKT_DIR <- file.path(BASE_DIR, "NKT")
dir.create(FIGURE_DIR, showWarnings = FALSE, recursive = TRUE)
setwd(FIGURE_DIR)

dpi <- 300

# ============ Part 1: Macrophage cytokine expression ============
cat("=== Part 1: Macrophage cytokine expression ===\n")

Macrophage.Integrated <- readRDS(file.path(MYELOID_DIR, "2-Macrophage.rds"))
cat("Loaded macrophage:", ncol(Macrophage.Integrated), "cells, clusters:",
    paste(levels(Idents(Macrophage.Integrated)), collapse=", "), "\n")

# Annotate macrophage groups based on our clustering
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
Macrophage.Integrated$macro_group <- Idents(Macrophage.Integrated)
Macrophage.Integrated <- subset(Macrophage.Integrated, subset = macro_group != 'Other')
macro_groups <- c('Group1', 'Group2', 'Group3', 'Group4')
Macrophage.Integrated$macro_group <- factor(Macrophage.Integrated$macro_group,
                                             levels = macro_groups)
Idents(Macrophage.Integrated) <- Macrophage.Integrated$macro_group

# Cytokines and chemokines
cytokines <- c('IL1B','IL6','IL10','TNF','IFNG',
               'CXCL1','CXCL2','CXCL8','CXCL9','CXCL10','CXCL11','CXCL16',
               'CCL2','CCL3','CCL4','CCL7','CCL8','CCL13','CCL18',
               'CX3CL1','XCL1','XCL2')
existing_cyto <- cytokines[cytokines %in% rownames(Macrophage.Integrated)]
cat("Plotting", length(existing_cyto), "/", length(cytokines), "cytokines\n")

# Cytokine dotplot per sample
DefaultAssay(Macrophage.Integrated) <- "RNA"
Idents(Macrophage.Integrated) <- 'sample_new'
samples_order <- c('HC1','HC2','HC3','HC4','O1','O2','O3','S1','C1','C2','C3','C4','C5')
Macrophage.Integrated@meta.data$sample_new <- factor(
  Macrophage.Integrated@meta.data$sample_new, levels = samples_order)

pdf(file = "5-heatmap-cytokine_sample.pdf", width = 11, height = 4)
pp <- DotPlot(Macrophage.Integrated, features = rev(existing_cyto),
              cols = c('white', '#F8766D'), dot.scale = 5) + RotatedAxis()
pp <- pp + theme(axis.text.x = element_text(size = 9), axis.text.y = element_text(size = 9)) +
  labs(x = '', y = '') +
  guides(color = guide_colorbar(title = 'Scale expression'),
         size = guide_legend(title = 'Percent expressed'))
print(pp)
dev.off()
cat("  -> 5-heatmap-cytokine_sample.pdf\n")

# Cytokine expression by macrophage group
Idents(Macrophage.Integrated) <- Macrophage.Integrated$macro_group
pdf(file = "5-heatmap-cytokine_macro_group.pdf", width = 8, height = 4)
pp <- DotPlot(Macrophage.Integrated, features = rev(existing_cyto),
              cols = c('white', '#F8766D'), dot.scale = 5) + RotatedAxis()
pp <- pp + theme(axis.text.x = element_text(size = 9), axis.text.y = element_text(size = 10)) +
  labs(x = '', y = '') +
  guides(color = guide_colorbar(title = 'Scale expression'),
         size = guide_legend(title = 'Percent expressed'))
print(pp)
dev.off()
cat("  -> 5-heatmap-cytokine_macro_group.pdf\n")

# Select cytokines for violin plots
cytokines_select <- c('IL1B','IL6','TNF','CXCL8','CXCL9','CXCL10',
                      'CCL2','CCL3','CCL4','CCL7','CCL8','CCL18')
existing_select <- cytokines_select[cytokines_select %in% rownames(Macrophage.Integrated)]

Idents(Macrophage.Integrated) <- 'group'
Macrophage_violin <- subset(Macrophage.Integrated, idents = c('O', 'S/C'))

pdf(file = "5-violin_cytokine.pdf", width = 12, height = 6)
vp_list <- VlnPlot(Macrophage_violin, features = existing_select,
                    pt.size = 0, combine = FALSE, ncol = 1)
vp_list <- lapply(vp_list, function(p) p + labs(x = '') +
  theme(axis.text = element_text(size = 8),
        axis.text.x = element_text(size = 8, angle = 0, hjust = 0.5),
        axis.title = element_text(size = 8),
        plot.title = element_text(size = 8, face = 'italic')))
print(patchwork::wrap_plots(vp_list, ncol = 6))
dev.off()
cat("  -> 5-violin_cytokine.pdf\n")

# Differential expression of cytokines between groups
Idents(Macrophage.Integrated) <- 'group'
cat("Differential cytokine expression (Wilcoxon)...\n")

cytokines_all <- unique(c(cytokines,
  'CXCL3','CXCL5','CXCL6','CXCL12','CXCL13','CXCL14','CXCL17',
  'CCL1','CCL5','CCL3L1','CCL4L1','CCL14','CCL15','CCL19','CCL20',
  'CCL22','CCL23','CCL24','CCL25','CCL26','CCL27','CCL28',
  'IL4','IL5','IL12B','IFNA1'))
cytokines_existing <- cytokines_all[cytokines_all %in% rownames(Macrophage.Integrated)]

deg1 <- FindMarkers(Macrophage.Integrated, ident.1 = 'S/C', ident.2 = 'O',
                     features = cytokines_existing, logfc.threshold = 0,
                     min.pct = 0, only.pos = FALSE)
deg1$comparison <- 'Severe_vs_Moderate'
deg1$gene <- rownames(deg1)
write.table(deg1, file = '5-cytokine-deg-severe-mild.txt', quote = FALSE, sep = '\t')

deg2 <- FindMarkers(Macrophage.Integrated, ident.1 = 'O', ident.2 = 'HC',
                     features = cytokines_existing, logfc.threshold = 0,
                     min.pct = 0, only.pos = FALSE)
deg2$comparison <- 'Moderate_vs_HC'
deg2$gene <- rownames(deg2)
write.table(deg2, file = '5-cytokine-deg-mild-hc.txt', quote = FALSE, sep = '\t')

deg3 <- FindMarkers(Macrophage.Integrated, ident.1 = 'S/C', ident.2 = 'HC',
                     features = cytokines_existing, logfc.threshold = 0,
                     min.pct = 0, only.pos = FALSE)
deg3$comparison <- 'Severe_vs_HC'
deg3$gene <- rownames(deg3)
write.table(deg3, file = '5-cytokine-deg-severe-hc.txt', quote = FALSE, sep = '\t')
cat("  -> 5-cytokine-deg-*.txt (3 files)\n")

# Free macrophage memory
rm(Macrophage.Integrated, Macrophage_violin)
gc()

# ============ Part 2: Chemokine receptor across cell types ============
cat("\n=== Part 2: Chemokine receptors across cell types ===\n")

# Load NKT annotated data first
NKT.clean <- readRDS(file.path(NKT_DIR, "6-NKT-annotated.rds"))
if ("celltype" %in% colnames(NKT.clean@meta.data)) {
  Idents(NKT.clean) <- NKT.clean$celltype
}
nkt_types <- unique(Idents(NKT.clean))
nkt_types <- nkt_types[nkt_types != 'Uncertain']
NKT.clean <- subset(NKT.clean, idents = nkt_types)
cat("NKT types:", paste(nkt_types, collapse=", "), "\n")

# Compute NKT averages now
DefaultAssay(NKT.clean) <- "RNA"
nkt_avg <- AverageExpression(NKT.clean, assays = "RNA", group.by = "celltype")$RNA
rm(NKT.clean)
gc()

# Load nCoV, extract other cell types, free it
nCoV.integrated <- readRDS(file.path(BASE_DIR, "nCoV.rds"))
nCoV.integrated <- RenameIdents(nCoV.integrated,
  '13'='Epithelial','16'='Epithelial','25'='Epithelial','28'='Epithelial','31'='Epithelial',
  '0'='Macrophages','1'='Macrophages','2'='Macrophages','3'='Macrophages',
  '4'='Macrophages','5'='Macrophages','7'='Macrophages','8'='Macrophages',
  '10'='Macrophages','11'='Macrophages','12'='Macrophages','18'='Macrophages',
  '21'='Macrophages','22'='Macrophages','23'='Macrophages','26'='Macrophages',
  '30'='Mast','6'='T','9'='T','14'='T','17'='NK','19'='Plasma',
  '27'='B','15'='Neutrophil','20'='mDC','29'='pDC','24'='Doublets')
nCoV.integrated$celltype_broad <- Idents(nCoV.integrated)
nCoV.integrated <- subset(nCoV.integrated, subset = celltype_broad != 'Doublets')
other_cells <- subset(nCoV.integrated, idents = c('Neutrophil', 'B', 'Plasma', 'pDC'))
DefaultAssay(other_cells) <- "RNA"
other_avg <- AverageExpression(other_cells, assays = "RNA", group.by = "celltype_broad")$RNA
rm(nCoV.integrated, other_cells)
gc()

# Reload macrophage for averages
Macrophage.Integrated <- readRDS(file.path(MYELOID_DIR, "2-Macrophage.rds"))
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
Macrophage.Integrated$celltype <- Idents(Macrophage.Integrated)
Macrophage.Integrated <- subset(Macrophage.Integrated, subset = celltype != 'Other')
DefaultAssay(Macrophage.Integrated) <- "RNA"
macro_avg <- AverageExpression(Macrophage.Integrated, assays = "RNA", group.by = "celltype")$RNA
rm(Macrophage.Integrated)
gc()

# Combine averages
chemokine_receptors <- c('CCR1','CCR2','CCR4','CCR5','CCR6','CCR7',
                          'CXCR1','CXCR2','CXCR3','CXCR4','CXCR5','CXCR6',
                          'CX3CR1','XCR1')
existing_receptors <- chemokine_receptors[chemokine_receptors %in% rownames(nkt_avg)]

common_genes <- intersect(intersect(rownames(macro_avg), rownames(nkt_avg)),
                          rownames(other_avg))
common_genes <- intersect(common_genes, existing_receptors)

combined_avg <- cbind(macro_avg[common_genes, , drop=FALSE],
                      nkt_avg[common_genes, , drop=FALSE],
                      other_avg[common_genes, , drop=FALSE])
cat("Combined matrix:", nrow(combined_avg), "receptors x", ncol(combined_avg), "cell types\n")

write.table(combined_avg, file = '5-chemokine-receptor-average.txt',
            quote = FALSE, sep = '\t', row.names = TRUE)
cat("  -> 5-chemokine-receptor-average.txt\n")

cat("\n=== Cytokine analysis complete ===\n")
cat("Files saved to:", FIGURE_DIR, "\n")
