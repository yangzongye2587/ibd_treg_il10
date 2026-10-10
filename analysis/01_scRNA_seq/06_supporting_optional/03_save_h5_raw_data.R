setwd('data')
getwd()
library(scKidney)
IBD_T_CELL_remove_dcs2 <- readRDS("data/IBD_T_CELL_remove_dcs2.rds")
DotPlot(IBD_T_CELL_remove_dcs2,features = c('CD4','CD8A','CD3E','TRAC','ITGA1'),group.by = 'T_subtype')
marker_genes_list <- list(
  "General T cell markers" = c('PTPRC','JCHAIN',"CD3D", "CD3E", "CD3G", "CD8A", "CD8B", "CD4",'TRAC'),
  "ILC3" = c("IL7R","RORC","GATA3","IL2RG","ID2","TOX","TCF7"),
  "MAIT" = c("SLC4A10", "TRAV1-2"),
  "NK_CD16" = c("FCGR3A", "GNLY", "GZMB", "NKG7"),
  "NK_CD56bright" = c("XCL1", "IL2RB", "NCAM1", "NCR1"),
  "T follicular helper markers" = c("CXCR5", "BCL6", "PDCD1", "CD40LG"),
  "Naive T cell markers" = c("CCR7", "SELL",'IL17A','IFNG'),
  "Tregs" = c("FOXP3", "CTLA4", "IL10"),
  "Trm/em_CD8" = c("GZMK", "CRTAM", "EOMES"),
  "Tissue resident memory T markers" = c("ITGAE", "ITGA1", "SPRY1",'CD69','RUNX3','CXCR6'),
  "gdT markers" = c("KLRC2", "TRDC", "TRGC1", "KIR2DL4"),
  "proli markers" = c("MKI67", "TOP2A", "PCNA"),
  'extro' = c('TIGIT', 'ICOS','LAG3'),
  'NKT' = c('TBX21','TNF')
)

DotPlot(
  object   = IBD_T_CELL_remove_dcs2,
  group.by = "T_subtype",
  features = marker_genes_list
) + RotatedAxis()
# 假设您的 Seurat 对象名为 seurat_object

# 指定输出的 H5ad 文件路径
library(scKidney)
output_h5ad_path <- "data\\treg_sub.h5ad"

major_annotation_ibd_subtype$cell_subtype.x <- major_annotation_ibd_subtype$cell_subtype
head(major_annotation_ibd_subtype@meta.data)
kidneyH5(TREG_SUBSET, output_h5ad_path)
CD4_IBD_T <- subset(IBD_T_CELL_remove_dcs2,T_subtype %in% c('Naive_T','TFH','TH17','Treg')&disease=='UC')
###v5

library(srtdisk)
library(Seurat)

obj <- IBD_T_CELL_remove_dcs2
DefaultAssay(obj) <- "RNA"

# 保存为 h5Seurat
SaveH5Seurat(obj, filename = "IBD_T_CELL_remove_dcs2.h5Seurat", overwrite = TRUE)

# 转 h5ad（会把 counts 放进 raw，data 放进 X，取决于实现/参数）
Convert("IBD_T_CELL_remove_dcs2.h5Seurat", dest = "h5ad", overwrite = TRUE)


####v4
library(Seurat)
library(SeuratDisk)

obj <- major_annotation_ibd_subtype 
DefaultAssay(obj) <- "RNA"

# 1) 从 v5 的 layer 里取出矩阵
cnt  <- GetAssayData(obj, assay = "RNA", layer = "counts")
dat  <- GetAssayData(obj, assay = "RNA", layer = "data")   # log-normalized

# 2) 重新创建一个“旧 Assay”
obj[["RNA"]] <- CreateAssayObject(counts = cnt)
obj[["RNA"]]@data <- dat

# （可选）减肥，避免写入巨大 scale.data
obj[["RNA"]]@scale.data <- matrix(, 0, 0)

# 3) 再走 SeuratDisk 的两步
SaveH5Seurat(obj, filename = "IBD_T_CELL_remove_dcs2.h5Seurat", overwrite = TRUE)
Convert("IBD_T_CELL_remove_dcs2.h5Seurat", dest = "h5ad", overwrite = TRUE)

