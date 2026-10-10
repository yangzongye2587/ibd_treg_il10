meta <- GSE214695@meta.data
head(meta)
meta$gsm <- sub("^(GSM[0-9]+).*", "\\1", rownames(meta))
meta$cell_name <- rownames(meta)
head(meta$gsm)
data$gsm <- rownames(data)
merge_meta <- meta %>% left_join(data,by = 'gsm')


head(merge_meta)
colnames(merge_meta)
merge_meta1 <- merge_meta %>% dplyr::select(c("orig.ident","nCount_RNA",
                                              "nFeature_RNA","group","gsm","cell_name","gender:ch1","tissue:ch1"))
colnames(merge_meta1)[7:8] <- c('sex','tissue')
merge_meta1$age <- NA
GSE214695@meta.data <- merge_meta1
gse <- getGEO("GSE231993", GSEMatrix = TRUE)
data <- gse[["GSE231993_series_matrix.txt.gz"]]@phenoData@data
meta <- merge_GSE231993@meta.data
meta$gsm <- sub("^(GSM[0-9]+).*", "\\1", rownames(meta))
meta$cell_name <- rownames(meta)
data$gsm <- data$geo_accession
merge_meta <- meta %>% left_join(data,by = 'gsm')
rownames(merge_meta) <- merge_meta$cell_name
merge_meta$group <- merge_meta$`sample group:ch1`
merge_meta$sex <- NA
merge_meta$tissue <- merge_meta$`tissue:ch1`
merge_meta1 <- merge_meta %>% dplyr::select(c("orig.ident","nCount_RNA",
                                              "nFeature_RNA","group","gsm","cell_name","sex","tissue"))
merge_meta1$age <- NA
merge_GSE231993@meta.data <- merge_meta1
colnames(META_DATA)
META_DATA <- GSE282122 @meta.data
META_DATA$orig.ident <- META_DATA$folder
META_DATA$sex <- META_DATA$`Sex:ch1`
META_DATA$tissue <- META_DATA$`site:ch1`
META_DATA$group <- META_DATA$`disease:ch1`
META_DATA$gsm <- META_DATA$geo_accession
META_DATA$cell_name <- rownames(META_DATA) 
META_DATA$age <- META_DATA$`age:ch1`
META_DATA1 <- META_DATA %>% dplyr::select(c("orig.ident","nCount_RNA",
                                              "nFeature_RNA","group","gsm","cell_name","sex","tissue",'age'))
GSE282122@meta.data <- META_DATA1
colnames(GSE282122@meta.data )
colnames(GSE214695@meta.data)
colnames(merge_GSE231993@meta.data)
# 1) 加来源列（后续做batch/分组统计很有用）
GSE282122$dataset <- "GSE282122"
GSE214695$dataset <- "GSE214695"
GSE231993$dataset <- "GSE231993"

# 2) 给细胞名加前缀，避免不同对象里有同名cell导致冲突
GSE282122 <- RenameCells(GSE282122, add.cell.id = "GSE282122")
GSE214695 <- RenameCells(GSE214695, add.cell.id = "GSE214695")
GSE231993 <- RenameCells(GSE231993, add.cell.id = "GSE231993")

# 3) 合并三个对象
merged_all <- merge(
  x = GSE282122,
  y = list(GSE214695, GSE231993),
  project = "Merged_3datasets"
)

# 4) 快速检查
dim(merged_all)
table(merged_all$dataset)
table(merged_all$dataset, merged_all$group)
colnames(META_DATA)
left_data <- META_DATA %>% select(c("orig.ident",'characteristics_ch1.7','characteristics_ch1.8','characteristics_ch1.14'))

merged_all_meta_data <- merged_all@meta.data

left_data_unique <- GSE231993_left_data %>% distinct(orig.ident, .keep_all = TRUE)
merged_all_meta_data$cell_name1 <- rownames(merged_all_meta_data) 
# 然后再进行合并

merged_all_meta_data <- merged_all_meta_data %>% left_join(left_data_unique, by = 'orig.ident')
rownames(merged_all_meta_data) <- merged_all_meta_data$cell_name1
merged_all@meta.data <- merged_all_meta_data
saveRDS(merged_all,file = 'uc_merge_data.rds')
saveRDS(GSE282122,file = 'GSE282122.rds')
saveRDS(GSE214695,file = 'GSE214695.rds')
saveRDS(merge_GSE231993,file = 'GSE231993.rds')
saveRDS(left_data,file = 'GSE231993_left_data.rds')
