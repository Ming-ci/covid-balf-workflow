.libPaths(c("E:/R/R_libs", .libPaths()))
library(Seurat)
nCoV.integrated <- readRDS("E:/Claude code/shengxin/covid_balf/nCoV.rds")
cat("=== NKT clusters (6,9,14,17) ===\n")
nkt <- subset(nCoV.integrated, idents = c('6','9','14','17'))
cat("NKT cells:", ncol(nkt), "\n")
tbl <- sort(table(nkt@meta.data$sample_new), decreasing = TRUE)
print(tbl)

cat("\n=== T cell clusters (0,1,2,3,5,6,9,10 from NKT=6,9) ===\n")
# T cells are from NKT: clusters 0,1,2,3,5,6,9,10
# But we need to first load NKT.rds to know these cluster IDs
# For now, check if 6-NKT.rds already exists
