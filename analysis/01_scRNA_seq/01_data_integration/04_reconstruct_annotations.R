IBD_MYLOID <- readRDS("data/IBD_MYLOID.rds")
IBD_T_CELL_remove_dcs2 <- readRDS("data/IBD_T_CELL_remove_dcs2.rds")
IMMUNE_MAJOR_IBD <- readRDS("data/IMMUNE_MAJOR_IBD.rds")
major_annotation <- readRDS("data/major_annotation.rds")
EPI_IBD <- readRDS("data/EPI_IBD.rds")
stromal_IBD <- readRDS("data/stromal_IBD.rds")

EPI <- EPI_IBD
stromal <- stromal_IBD
immune <- IMMUNE_MAJOR_IBD

DimPlot(major_annotation,group.by = 'RNA_snn_res.0.8',label = T)
major <- major_annotation@meta.data
IMMUNE <- immune@meta.data
Epi <- EPI@meta.data
Stromal <- stromal@meta.data
major$cell_subtype <- NA

IMMUNE$cell_subtype <- IMMUNE$celltype_immune
Epi$cell_subtype <- Epi$marjor_type
Stromal$cell_subtype <- Stromal$celltype_stromal
head(major)
a <- colnames(Epi)
a <- a['cell_subtype']
Epi <- Epi %>% dplyr::select('cell_subtype')
Stromal <- Stromal%>% select('cell_subtype')
IMMUNE <- IMMUNE%>% select('cell_subtype')
merge_data <- rbind(Epi,Stromal,IMMUNE)
merge_data$names <- rownames(merge_data)
major$names <- rownames(major)
merge_data <- major %>% left_join(merge_data,by = 'names')



head(merge_data)
# 首先，确保 merge_data 有正确的行名
rownames(merge_data) <- merge_data$names

major_annotation@meta.data <- merge_data
major_annotation$cell_subtype <- major_annotation$cell_subtype.y
DimPlot(major_annotation,group.by = 'cell_subtype')
saveRDS(major_annotation,file = 'major_annotation_ibd_subtype.rds')


library(tidyverse)
library(Seurat)
NEC_MAJORTYPE<- EDNRB_majortype2.0
NEC_MAJORTYPE1 <- NEC_MAJORTYPE
major <- NEC_MAJORTYPE@meta.data
IMMUNE <- immune@meta.data
Epi <- EPI@meta.data
Stromal <- stromal@meta.data
major$cell_subtype <- NA

IMMUNE$cell_subtype <- IMMUNE$celltype_immune
Epi$cell_subtype <- Epi$celltype_epi
Stromal$cell_subtype <- Stromal$celltype_Stromal
head(major)
a <- colnames(Epi)
a <- a[-21]
Epi <- Epi %>% dplyr::select(a)
Stromal <- Stromal%>% select(a)
IMMUNE <- IMMUNE%>% select(a)
merge_data <- rbind(Epi,Stromal,IMMUNE)
merge_data$names <- rownames(merge_data)
head(merge_data)
head(major)
NEC_MAJORTYPE$names <- rownames(major)

major_names <- data.frame(NEC_MAJORTYPE$names)
major_names1 <- major_names %>% left_join(merge_data)



head(major_names)
head(merge_data)
# 首先，确保 merge_data 有正确的行名
rownames(merge_data) <- merge_data$names

# 把 merge_data 中的注释向量提取出来
cell_subtype_vector <- merge_data$cell_subtype
names(cell_subtype_vector) <- merge_data$names

# 从 major_names 提取细胞名
cell_names <- rownames(major_names)

# 按顺序提取对应的注释（如果不存在则为 NA）
major_names$merge_sub.cell_subtype <- cell_subtype_vector[cell_names]
NEC_MAJORTYPE$cell_subtype <- NA
NEC_MAJORTYPE$cell_subtype <- major_names$merge_sub.cell_subtype


DotPlot(IBD_T_CELL_remove_dcs2,features = c('CD4','CD8A'),group.by = 'T_subtype')
DimPlot(NEC_MAJORTYPE,reduction = 'umap',group.by = 'cell_subtype',label = T,cols = c(cell_colors,clusterCols,cell_colors))