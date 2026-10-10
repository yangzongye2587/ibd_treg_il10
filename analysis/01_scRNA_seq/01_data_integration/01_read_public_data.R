setwd("data\\UC")  

library(tidyverse)
library(Seurat)
library(stringr)
data_dir <- 'GSE214695_RAW/' 
fs <- list.files(data_dir)

# 找到所有 matrix 文件（每个样本一个）
mtx_files <- fs[str_detect(fs, "_matrix\\.mtx\\.gz$")]

# 从文件名提取样本ID：去掉 "_matrix.mtx.gz"
sample_ids <- str_replace(mtx_files, "_matrix\\.mtx\\.gz$", "")

# 逐个样本读取，生成 Seurat 对象列表
objs <- lapply(sample_ids, function(sid){
  mtx <- file.path(data_dir, paste0(sid, "_matrix.mtx.gz"))
  feat <- file.path(data_dir, paste0(sid, "_features.tsv.gz"))
  bar  <- file.path(data_dir, paste0(sid, "_barcodes.tsv.gz"))
  
  counts <- ReadMtx(mtx = mtx, features = feat, cells = bar, feature.column = 2)
  
  obj <- CreateSeuratObject(counts = counts, project = sid, min.cells = 3, min.features = 200)
  obj$sample_id <- sid
  
  # 如果你想从sid里再解析分组（HC/UC/CD）
  obj$group <- str_extract(sid, "(HC|UC|CD)")
  obj
})

names(objs) <- sample_ids

# 合并成一个大对象（后续整合/注释用）
obj_all <- merge(objs[[1]], y = objs[-1], add.cell.ids = names(objs))
table(obj_all$group)
obj_all <- subset(obj_all,group %in% c('HC',    'UC'))

saveRDS(obj_all,file = 'GSE214695.rds')

library(tidyverse)
library(Seurat)
library(stringr)
data_dir <- 'GSE231993_RAW/' 
fs <- list.files(data_dir)

# 找到所有 matrix 文件（每个样本一个）
mtx_files <- fs[str_detect(fs, "-matrix\\.mtx\\.gz$")]

# 从文件名提取样本ID：去掉 "_matrix.mtx.gz"
sample_ids <- str_replace(mtx_files, "-matrix\\.mtx\\.gz$", "")

# 逐个样本读取，生成 Seurat 对象列表
objs1 <- lapply(sample_ids, function(sid){
  mtx <- file.path(data_dir, paste0(sid, "-matrix.mtx.gz"))
  feat <- file.path(data_dir, paste0(sid, "-features.tsv.gz"))
  bar  <- file.path(data_dir, paste0(sid, "-barcodes.tsv.gz"))
  
  counts <- ReadMtx(mtx = mtx, features = feat, cells = bar, feature.column = 2)
  
  obj <- CreateSeuratObject(counts = counts, project = sid, min.cells = 3, min.features = 200)
  obj$sample_id <- sid
  
  # 如果你想从sid里再解析分组（HC/UC/CD）
  obj$group <- str_extract(sid, "(HC|UC|CD)")
  obj
})

names(objs1) <- sample_ids

# 合并成一个大对象（后续整合/注释用）
obj_all1 <- merge(objs1[[1]], y = objs1[-1], add.cell.ids = names(objs1))
obj_all1_all <- merge(obj_all, y = obj_all1)
unique(obj_all1_all$orig.ident)

saveRDS(obj_all1,file = 'merge_GSE231993.rds')

library(Seurat)
library(SeuratDisk)
library(Matrix)

# Define the directory where the folders are located
data_dir <- "data/UC/GSE282122_raw_processed_data/raw_processed_data"

# List all subfolders
fs <- list.files(data_dir, full.names = TRUE)

# Initialize a list to store Seurat objects
seurat_objects <- list()

# Iterate through all folders
for (folder in fs) {
  # Read the H5 file for the current folder
  file_path <- file.path(folder, "raw_feature_bc_matrix.h5")
  count <- Read10X_h5(file_path)
  
  # Create a Seurat object
  seurat_obj <- CreateSeuratObject(counts = count,
                                   project = "responder",
                                   min.cells = 5,
                                   min.features = 300)
  
  # Add the folder name as a metadata column
  seurat_obj$folder <- basename(folder)  # The folder name is added to the metadata
  
  # Store the Seurat object in the list
  seurat_objects[[basename(folder)]] <- seurat_obj
}
rm(seurat_objects)
# Merge all Seurat objects into one
combined_seurat_object <- merge(seurat_objects[[1]], y = seurat_objects[2:length(seurat_objects)], add.cell.ids = names(seurat_objects))

# You can now explore the combined Seurat object
head(combined_seurat_object@meta.data)


library(GEOquery)
gse <- getGEO("GSE282122", GSEMatrix = TRUE)
# 获取表达数据和临床信息
expression_data <- exprs(gse[[1]])  # 基因表达矩阵
pheno_data <- pData(gse[[1]])       # 临床信息
table(pheno_data$`disease:ch1`)
pheno_data1 <- pheno_data %>% filter(`disease:ch1` %in% c('UC','Healthy'))
filter_vector <- pheno_data1$title
pheno_data1$title_raw <- pheno_data1$title
pheno_data1$title <- sub("-reup$", "", pheno_data1$title)

table(pheno_data1$title)
table(combined_seurat_object$folder)
h5_sub <- subset(combined_seurat_object,folder %in% c(unique(pheno_data1$title)))
h5_sub$title <- h5_sub$folder
meta <-h5_sub@meta.data
meta$cell_name <- rownames(meta)
metamerge <- meta %>% left_join(pheno_data1,by = 'title')
metamerge <- metamerge %>% column_to_rownames('cell_name')
h5_sub@meta.data <- metamerge
head(h5_sub@meta.data )
saveRDS(h5_sub,file = 'GSE282122.rds')
setwd("data\\UC")  
table(uc_merge_data$characteristics_ch1.7,uc_merge_data$characteristics_ch1.8)
table(is.na(uc_merge_data$characteristics_ch1.14))
# 如果是 character 向量
uc_merge_data$characteristics_ch1.8[is.na(uc_merge_data$characteristics_ch1.8)] <- "treatment: None"
uc_merge_data_disease <- subset(uc_merge_data, characteristics_ch1.8 != "treatment: Post")
table(uc_merge_data_disease$characteristics_ch1.7,uc_merge_data_disease$characteristics_ch1.8)
table(is.na(uc_merge_data_disease$characteristics_ch1.14))
uc_merge_data_disease$disease <- ifelse(uc_merge_data_disease$group %in% c('HC','Healthy'),
                                        'HC','UC')
uc_merge_data_disease$sex1 <- ifelse(uc_merge_data_disease$sex %in% c('F','Female'),
                                     'Female','Male')


table(uc_merge_data_disease$sex1)
saveRDS(uc_merge_data_disease,file = 'uc_merge_data_disease_filterover.rds')
