library(scKidney)

# 假设您的 Seurat 对象名为 seurat_object

# 指定输出的 H5ad 文件路径

output_h5ad_path <- "UC_treg_IBD_3group.h5ad"
setwd("data\\UC")  
# 运行转换函数
getwd()
sce.all.filt <- TREG_SUBSET
kidneyH5(sce.all.filt, output_h5ad_path)
sce.all.filt$cell_subtype<- sce.all.filt$T_subtype
DGE_results <- list()
for (tm in unique(sce.all.filt$cell_subtype)) {
  
  # 按细胞类型 subset
  nec_data_subset <- subset(sce.all.filt, subset = cell_subtype == tm)
  if (ncol(nec_data_subset) == 0) next
  
  # 设定分组为 age_disease
  Idents(nec_data_subset) <- "disease"
  

  
  # NEC vs Con 的差异分析（上下调都保留）
  dge <- FindMarkers(
    nec_data_subset,
    ident.1 = "UC",
    ident.2 = "HC",
    logfc.threshold = 0.25,   # 🔑 改小
    min.pct = 0.05,          # 🔑 改小
    min.diff.pct = 0.1,      # 🔑 改小
    test.use = "wilcox",
    only.pos = FALSE,
    assay = "RNA"
  )
  
  dge$celltype  <- tm
  
  # 标记方向（NEC 上调 or Con 上调）
  fc_col <- if ("avg_log2FC" %in% colnames(dge)) "avg_log2FC" else "avg_logFC"
  dge$direction <- ifelse(dge[[fc_col]] > 0, "NEC_p0_up", "Con_p0_up")
  dge <- dge %>% filter(p_val<0.05)
  dge$gene <- rownames(dge)
  DGE_results[[paste0(tm, "_UC_vs_Con")]] <- dge
}
