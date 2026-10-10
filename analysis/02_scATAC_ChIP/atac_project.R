options(BioC_mirror="https://mirrors.westlake.edu.cn/bioconductor")
dir <- 'data/scATAC'
setwd(dir)
samples=list.files( dir )
#BiocManager::install("sparseMatrixStats")
devtools::install_github("const-ae/sparseMatrixStats")
getwd()
BiocManager::install("EnsDb.Hsapiens.v75", force = TRUE)
untar("EnsDb.Hsapiens.v75_2.99.0.tar.gz",
      files = "EnsDb.Hsapiens.v75/inst/extdata/EnsDb.Hsapiens.v75.sqlite")
install.packages(  "data/scATAC\\EnsDb.Hsapiens.v75_2.99.0.tar.gz",
                        force = TRUE )




install.packages(  "data/reference\\EnsDb.Hsapiens.v75_2.99.0.tar.gz", 
                   repos = NULL,
                   type = "source"  )
library(ensembldb)
edb <- EnsDb("EnsDb.Hsapiens.v75/inst/extdata/EnsDb.Hsapiens.v75.sqlite")





dir <- 'data/scATAC'
setwd(dir)
samples=list.files( dir )
#library(EnsDb.Hsapiens.v75)  ##这个包用于人类的基因组注释文件
library(patchwork)
library(Signac)
library(Seurat)
library(ggplot2)
library(hdf5r)
library(spatstat.geom)
library(biovizBase)
library(Matrix)
library(readr)
library(tidyverse)
library(BiocParallel)
register(SerialParam())
library(chromVAR)
library(Signac)
library(Seurat)
library(chromVAR)
library(SummarizedExperiment)
library(BSgenome.Hsapiens.UCSC.hg38)
library(motifmatchr)
library(GenomeInfoDb)
library(Matrix)

metadata <- read_csv("scATAC_cell_metadata.csv")
counts <- readMM("scATAC_peak_matrix_tcell.mtx")

peaks <- read.table(
  "scATAC_peak_matrix_row_peaks_tcell.txt",
  stringsAsFactors = FALSE
)[,1]

cells <- read.table("data/cd_hu_data/scATAC_peak_matrix_col_cellIDs_tcell.txt",
                                                   quote="\"", comment.char="")[,1]

rownames(counts) <- peaks
colnames(counts) <- cells

chrom_assay <- CreateChromatinAssay(
  counts = counts,
  sep = c(":", "-"),
  genome = "hg38",
  min.cells = 10,
  min.features = 200
)

atac <- CreateSeuratObject(
  chrom_assay,
  assay = "peaks",
  project = "scATAC_Tcell"
)




# meta处理
meta_at <- atac@meta.data
meta_at$cell_id <- rownames(meta_at)
metadata$cell_id <- metadata$...1
meta_merge <- left_join(meta_at,metadata,by = 'cell_id')
meta_merge <- meta_merge %>% column_to_rownames('...1')
table(meta_merge$Patient_ID)
table(meta_merge$biopsy_location)
table(meta_merge$Seq_Library_ID)
meta_merge$group <- ifelse(meta_merge$Patient_ID %in% c('C1','C2','C3','C4','C5','C6'),
                           'Healthy','UC')
meta_merge$group_location <- paste0(meta_merge$group,'_',meta_merge$biopsy_location)
atac@meta.data <- meta_merge

atac <- RunTFIDF(atac)
atac <- FindTopFeatures(atac, min.cutoff = "q0")
atac <- RunSVD(atac)
atac <- RunUMAP(atac, reduction = "lsi", dims = 2:30)
atac <- FindNeighbors(atac, reduction = "lsi", dims = 2:30)
atac <- FindClusters(atac, verbose = FALSE)

# 读取 seqinfo
setMethod(
  "seqinfo",
  signature(x = "ChromatinAssay"),
  function(x) {
    GenomeInfoDb::seqinfo(x@ranges)
  }
)

# 写入 seqinfo
setReplaceMethod(
  "seqinfo",
  signature(x = "ChromatinAssay"),
  function(x, value) {
    GenomeInfoDb::seqinfo(x@ranges) <- value
    x
  }
)

# 读取 genome
setMethod(
  "genome",
  signature(x = "ChromatinAssay"),
  function(x) {
    GenomeInfoDb::genome(x@ranges)
  }
)

# 写入 genome
setReplaceMethod(
  "genome",
  signature(x = "ChromatinAssay"),
  function(x, value) {
    GenomeInfoDb::genome(x@ranges) <- value
    x
  }
)
GenomeInfoDb::genome(atac[["peaks"]])
genome(atac[["peaks"]]) <- "hg38"
GenomeInfoDb::genome(atac[["peaks"]])
annotations <- Annotation(atac) 
saveRDS(atac,file = 'atac3.rds')



Annotation(atac) <- annotations

cellRatioPlot(object = atac,
              sample.name = "group",
              celltype.name = "cluster_label3",fill.col = cell_colors)
DimPlot(atac,group.by =  "cluster_label3",cols =cell_colors)








options(Seurat.object.assay.version = "v5")
library(JASPAR2020)
library(TFBSTools)
library(BSgenome.Hsapiens.UCSC.hg38)

pfm <- getMatrixSet(
  JASPAR2020,
  opts = list(
    collection = "CORE",
    tax_group = "vertebrates",
    all_versions = FALSE
  )
)

atac <- AddMotifs(
  object = atac,
  genome = BSgenome.Hsapiens.UCSC.hg38,
  pfm = pfm
)
library(BiocParallel)
register(SerialParam())
library(chromVAR)
library(Signac)
library(Seurat)
library(chromVAR)
library(SummarizedExperiment)
library(BSgenome.Hsapiens.UCSC.hg38)
library(motifmatchr)
library(GenomeInfoDb)
library(Matrix)

# 1. 取 peak count matrix
peak_counts <- GetAssayData(
  object = atac,
  assay = "peaks",
  layer = "counts"
)


# 2. 从 peak 名字生成 GRanges
# 你的 peak 名字应该是 chr1:100-200 这种格式
peak_gr <- StringToGRanges(
  regions = rownames(peak_counts),
  sep = c(":", "-")
)

seqlevelsStyle(peak_gr) <- "UCSC"
genome(peak_gr) <- "hg38"

# 3. 构建 SummarizedExperiment
se <- SummarizedExperiment(
  assays = list(counts = peak_counts),
  rowRanges = peak_gr,
  colData = atac@meta.data[colnames(peak_counts), , drop = FALSE]
)

# 4. 加 GC bias
se <- addGCBias(
  object = se,
  genome = BSgenome.Hsapiens.UCSC.hg38
)

# 5. 用 JASPAR motif 匹配 peaks
# 如果你前面已经有 pfm，就不用重新 getMatrixSet
motif_ix <- motifmatchr::matchMotifs(
  pwms = pfm,
  subject = se,
  genome = BSgenome.Hsapiens.UCSC.hg38
)

# 6. 计算 chromVAR deviation score
dev <- computeDeviations(
  object = se,
  annotations = motif_ix
)

# 7. 提取 z-score matrix
dev_mat <- deviationScores(dev)

# 8. 放回 Seurat 对象
rownames(dev_mat) <- paste0("chromvar_", rownames(dev_mat))

atac[["chromvar"]] <- CreateAssayObject(data = dev_mat)
grep("FOXP3", rownames(atac[["chromvar"]]), value = TRUE, ignore.case = TRUE)
AAA <- rownames(atac[["chromvar"]])
library(Signac)
library(chromVAR)
# 看看 BATF 对应哪个 motif ID
ConvertMotifID(
  object = atac,
  name = "BATF"
)
batf_id <- ConvertMotifID(atac, name = "BATF")
batf_id

FeaturePlot(
  object = atac,
  features = batf_id,
  reduction = "umap"
)
atac$cell_group <- paste0(atac$cluster_label3,'_',atac$group)
VlnPlot(
  object = atac,
  features = batf_id,
  group.by = "cell_group",
  pt.size = 0
)
ids <- rownames(atac[["chromvar"]])

motif_names <- ConvertMotifID(
  object = atac,
  id = ids
)

head(data.frame(id = ids, name = motif_names))



atac$cell_group <- paste0(atac$cluster_label3,'_',atac$group)
batf_id <- ConvertMotifID(atac, name = "BATF")
batf_id

VlnPlot(atac, features = batf_id, group.by = "cell_group", pt.size = 0)

df <- FetchData(atac, vars = c(batf_id, "cell_group"))

aggregate(df[[batf_id]], by = list(cell_group = df$cell_group), FUN = mean)
library(dplyr)
colnames(df)
library(dplyr)
library(dplyr)
library(tidyr)

batf_col <- "chromvar_MA1634.1"

library(dplyr)

batf_col <- "chromvar_MA1634.1"
library(Seurat)
library(Signac)
library(GenomeInfoDb)
library(GenomicRanges)
#obj_sub <- subset(atac, subset = biopsy_location %in% c('AC','PP','TI'))

df <- FetchData(
  atac,
  vars = c(batf_col, "cluster_label3", "group", "biopsy_location")
)
library(dplyr)
library(ggplot2)

batf_col <- "chromvar_MA1634.1"

# 只保留 UC 和 Healthy
df_plot <- df %>%
  filter(group %in% c("Healthy", "UC"))

# 重新计算统计
batf_test <- df_plot %>%
  group_by(cluster_label3) %>%
  summarise(
    n_UC = sum(group == "UC"),
    n_Healthy = sum(group == "Healthy"),
    mean_UC = mean(.data[[batf_col]][group == "UC"], na.rm = TRUE),
    mean_Healthy = mean(.data[[batf_col]][group == "Healthy"], na.rm = TRUE),
    median_UC = median(.data[[batf_col]][group == "UC"], na.rm = TRUE),
    median_Healthy = median(.data[[batf_col]][group == "Healthy"], na.rm = TRUE),
    p_value = tryCatch(
      wilcox.test(.data[[batf_col]] ~ group)$p.value,
      error = function(e) NA_real_
    ),
    .groups = "drop"
  ) %>%
  mutate(
    diff = mean_UC - mean_Healthy,
    p_adj = p.adjust(p_value, method = "BH")
  ) %>%
  arrange(p_adj)

# 星号标签
p_labels <- batf_test %>%
  mutate(
    label = case_when(
      p_adj < 0.001 ~ "***",
      p_adj < 0.01  ~ "**",
      p_adj < 0.05  ~ "*",
      TRUE ~ "ns"
    )
  )

# 每个分面的星号位置
y_pos <- df_plot %>%
  group_by(cluster_label3) %>%
  summarise(y = max(.data[[batf_col]], na.rm = TRUE) * 1.05)

p_labels <- left_join(p_labels, y_pos, by = "cluster_label3")
p_labels$x <- 1.5

# 固定顺序
df_plot$group <- factor(df_plot$group, levels = c("Healthy", "UC"))

# 作图
ggplot(df_plot, aes(x = group, y = .data[[batf_col]], fill = group)) +
  geom_violin(trim = FALSE) +
  geom_boxplot(width = 0.15, outlier.size = 0.2) +
  facet_wrap(~ cluster_label3, scales = "free_y") +
  geom_text(
    data = p_labels,
    aes(x = x, y = y, label = label),
    inherit.aes = FALSE,
    size = 5
  ) +
  xlab("Group") +
  ylab("BATF chromVAR score")
cluster_order <- c("Treg", setdiff(unique(df_plot$cluster_label3), "Treg"))
df_plot$cluster_label3 <- factor(df_plot$cluster_label3, levels = cluster_order)
p_labels$cluster_label3 <- factor(p_labels$cluster_label3, levels = cluster_order)
p_labels2 <- p_labels %>% filter(label != "ns")

ggplot(df_plot, aes(x = group, y = .data[[batf_col]], fill = group)) +
  geom_violin(trim = FALSE) +
  geom_boxplot(width = 0.15, outlier.size = 0.2) +
  facet_wrap(~ cluster_label3, scales = "free_y") +
  geom_text(
    data = p_labels2,
    aes(x = x, y = y, label = label),
    inherit.aes = FALSE,
    size = 5
  )+
  xlab("Group") +
  ylab("BATF chromVAR score")

FeaturePlot(
  atac,
  features = batf_col,
  reduction = "umap"
)
cells_treg <- rownames(atac@meta.data)[atac$cluster_label3 == "Treg"]

DimPlot(atac, group.by = "cluster_label3")
library(dplyr)
library(ggplot2)

df_summary <- df %>%
  group_by(cluster_label3, group) %>%
  summarise(
    mean_batf = mean(.data[[batf_col]], na.rm = TRUE),
    .groups = "drop"
  )

ggplot(df_summary, aes(x = cluster_label3, y = mean_batf, fill = group)) +
  geom_col(position = "dodge") +
  coord_flip() +
  xlab("Cell type") +
  ylab("Mean BATF chromVAR score")
library(dplyr)
library(tidyr)
library(pheatmap)

mat_batf <- df %>%
  group_by(cluster_label3, group) %>%
  summarise(mean_batf = mean(.data[[batf_col]], na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = group, values_from = mean_batf) %>%
  as.data.frame()

rownames(mat_batf) <- mat_batf$cluster_label3
mat_batf$cluster_label3 <- NULL

pheatmap(as.matrix(mat_batf), scale = "none")




cells_treg <- rownames(atac@meta.data)[
  atac$cluster_label3 == "Treg" & atac$group %in% c("UC", "Healthy")
]
library(Matrix)

chromvar_mat <- GetAssayData(atac, assay = "chromvar", slot = "data")
chromvar_treg <- chromvar_mat[, cells_treg, drop = FALSE]
dim(chromvar_treg)
meta_treg <- atac@meta.data[cells_treg, c("group"), drop = FALSE]
table(meta_treg$group)
motif_test <- apply(chromvar_treg, 1, function(x) {
  uc <- x[meta_treg$group == "UC"]
  hc <- x[meta_treg$group == "Healthy"]
  
  c(
    mean_UC = mean(uc, na.rm = TRUE),
    mean_Healthy = mean(hc, na.rm = TRUE),
    diff = mean(uc, na.rm = TRUE) - mean(hc, na.rm = TRUE),
    p_value = tryCatch(wilcox.test(uc, hc)$p.value, error = function(e) NA_real_)
  )
})

motif_test <- as.data.frame(t(motif_test))
motif_test$motif_id <- rownames(motif_test)
motif_test$p_adj <- p.adjust(motif_test$p_value, method = "BH")
motif_test$motif_name <- ConvertMotifID(atac, id = motif_test$motif_id)
head(motif_test)
library(ggplot2)
library(dplyr)

motif_test$log10p <- -log10(motif_test$p_adj + 1e-300)
motif_test$sig <- ifelse(motif_test$p_adj < 0.05 & abs(motif_test$diff) > 0.2, "sig", "ns")

ggplot(motif_test, aes(x = diff, y = log10p, color = sig)) +
  geom_point(size = 1) +
  xlab("Mean chromVAR score difference (UC - Healthy)") +
  ylab("-log10 adjusted P")
library(ggrepel)

batf_row <- motif_test[grep("BATF", motif_test$motif_name, ignore.case = TRUE), ]

ggplot(motif_test, aes(x = diff, y = log10p, color = sig)) +
  geom_point(size = 1) +
  geom_point(data = batf_row, color = "red", size = 2) +
  geom_text_repel(data = batf_row, aes(label = motif_name), color = "red") +
  xlab("Mean chromVAR score difference (UC - Healthy)") +
  ylab("-log10 adjusted P")
top_motifs <- motif_test %>%
  filter(p_adj < 0.05) %>%
  arrange(desc(diff)) %>%
  head(30)

ggplot(top_motifs, aes(x = reorder(motif_name, diff), y = diff)) +
  geom_col(fill = '#4C78A8') +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  coord_flip() +
  xlab("Motif") +
  ylab("Mean chromVAR score difference (UC - Healthy)")

