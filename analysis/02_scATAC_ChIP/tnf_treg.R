############################################################
## BATF regulates TNFR2/TNFRSF1B in Treg - scATAC evidence
## Based on ArchR project: projHeme2
############################################################

suppressPackageStartupMessages({
  library(ArchR)
  library(dplyr)
  library(ggplot2)
  library(SummarizedExperiment)
  library(Matrix)
  library(GenomicRanges)
  library(GenomeInfoDb)
  library(rtracklayer)
  library(BSgenome.Hsapiens.UCSC.hg38)
})

addArchRThreads(threads = 1)

############################################################
## 0. 输入参数
############################################################
library(ArchR)
library(dplyr)

projHeme2 <- loadArchRProject(
  path = "data/IBD_ATAC/arrow/projHeme2_RNAIntegrated",
  force = TRUE
)

meta <- getCellColData(ArchRProj = projHeme2)
colnames(meta)
# 2. 提取整合后的 metadata
meta_df <- as.data.frame(getCellColData(ArchRProj = projHeme2))

# 加入细胞名
meta_df$cellNames <- rownames(meta_df)

# 3. 查看预测细胞类型 / 分群 / 置信度
head(meta_df[, c(
  "cellNames",
  "predictedCell_Un",
  "predictedGroup_Un",
  "predictedScore_Un"
)])
library(ArchR)
library(dplyr)

# 提取 metadata
meta_df <- as.data.frame(getCellColData(ArchRProj = projHeme2))
meta_df$cellNames <- rownames(meta_df)

# 根据 Sample 新建 Group 分组
meta_df$Group <- NA

meta_df$Group[grepl("Control", meta_df$Sample, ignore.case = TRUE)] <- "CON"
meta_df$Group[grepl("-PP", meta_df$Sample, ignore.case = TRUE)] <- "PP"

# 检查分组是否正确
table(meta_df$Sample, meta_df$Group, useNA = "ifany")
table(meta_df$Group, useNA = "ifany")
projHeme2 <- addCellColData(
  ArchRProj = projHeme2,
  data = meta_df$Group,
  cells = meta_df$cellNames,
  name = "Group",
  force = TRUE
)
meta_df <- as.data.frame(getCellColData(ArchRProj = projHeme2))
meta_df$cellNames <- rownames(meta_df)

table(meta_df$Sample, meta_df$Group)
treg_subtypes <- c(
  "actived_TREG",
  "NAIVE_LIKE_TREG",
  "TNF_TREG"
)
treg_meta <- meta_df %>%
  filter(predictedGroup_Un %in% treg_subtypes)

table(treg_meta$predictedGroup_Un)
table(treg_meta$Group, treg_meta$predictedGroup_Un)
treg_count <- as.data.frame(
  table(
    Group = treg_meta$Group,
    Treg_subtype = treg_meta$predictedGroup_Un
  )
)

treg_count

meta_df <- as.data.frame(getCellColData(projHeme2))
meta_df$cellNames <- rownames(meta_df)

# 根据 Sample 新建 CON / PP
meta_df$Group <- ifelse(
  grepl("Control", meta_df$Sample, ignore.case = TRUE),
  "CON",
  ifelse(grepl("-PP", meta_df$Sample, ignore.case = TRUE), "PP", NA)
)

table(meta_df$Sample, meta_df$Group, useNA = "ifany")

# Treg 亚群
treg_subtypes <- c(
  "actived_TREG",
  "NAIVE_LIKE_TREG",
  "TNF_TREG"
)

# 新建 Treg_Group

meta <- getCellColData(ArchRProj = projHeme2)
colnames(meta)
# 2. 提取整合后的 metadata
meta_df <- as.data.frame(getCellColData(ArchRProj = projHeme2))

# 加入细胞名
meta_df$cellNames <- rownames(meta_df)



meta_df$Treg_Group <- NA

sel <- meta_df$predictedGroup_Un %in% treg_subtypes & !is.na(meta_df$Group)

meta_df$Treg_Group[sel] <- paste(
  meta_df$predictedGroup_Un[sel],
  meta_df$Group[sel],
  sep = "_"
)

table(meta_df$Treg_Group, useNA = "ifany")
projHeme2 <- addCellColData(
  ArchRProj = projHeme2,
  data = meta_df$Group,
  cells = meta_df$cellNames,
  name = "Group",
  force = TRUE
)

projHeme2 <- addCellColData(
  ArchRProj = projHeme2,
  data = meta_df$Treg_Group,
  cells = meta_df$cellNames,
  name = "Treg_Group",
  force = TRUE
)

treg_cells <- meta_df$cellNames[!is.na(meta_df$Treg_Group)]


############################################################
## 0. 输入参数
############################################################



outDir <- "data/IBD_ATAC/BATF_TNFRSF1B_Treg_analysis_treg_sub"
dir.create(outDir, recursive = TRUE, showWarnings = FALSE)

# 这里改成你的 BATF ChIP-seq peak 文件路径
batf_chip_bed <- "data/IBD_ATAC/BATF_ChIPseq_peaks.bed"

target_gene <- "TNFRSF1B"   # TNFR2 的基因名
tf_gene <- "BATF"

treg_subtypes <- c(
  "actived_TREG",
  "NAIVE_LIKE_TREG",
  "TNF_TREG"
)

############################################################
## 1. 整理 metadata：新建 CON / PP 分组和 Treg_Group
############################################################

meta_df <- as.data.frame(getCellColData(projHeme2))
meta_df$cellNames <- rownames(meta_df)

meta_df$Group <- ifelse(
  grepl("Control", meta_df$Sample, ignore.case = TRUE),
  "CON",
  ifelse(grepl("-PP", meta_df$Sample, ignore.case = TRUE), "PP", NA)
)

meta_df$Treg_Group <- NA

sel_treg <- meta_df$predictedGroup_Un %in% treg_subtypes & !is.na(meta_df$Group)

meta_df$Treg_Group[sel_treg] <- paste(
  meta_df$predictedGroup_Un[sel_treg],
  meta_df$Group[sel_treg],
  sep = "_"
)

cat("\n===== Sample vs Group =====\n")
print(table(meta_df$Sample, meta_df$Group, useNA = "ifany"))

cat("\n===== Treg_Group =====\n")
print(table(meta_df$Treg_Group, useNA = "ifany"))

write.csv(
  meta_df,
  file = file.path(outDir, "01_all_cell_metadata_with_Group_TregGroup.csv"),
  row.names = FALSE
)

projHeme2 <- addCellColData(
  ArchRProj = projHeme2,
  data = meta_df$Group,
  cells = meta_df$cellNames,
  name = "Group",
  force = TRUE
)

projHeme2 <- addCellColData(
  ArchRProj = projHeme2,
  data = meta_df$Treg_Group,
  cells = meta_df$cellNames,
  name = "Treg_Group",
  force = TRUE
)

############################################################
## 2. 提取 Treg 子项目
############################################################

library(ArchR)
library(dplyr)
library(ggplot2)
library(SummarizedExperiment)

# 提取 metadata
meta_df <- as.data.frame(getCellColData(
  projHeme2,
  select = c("Sample", "Group", "predictedGroup_Un")
))

meta_df$cell <- rownames(meta_df)

# 只取 Treg，保留 CON 和 PP 两组
treg_cells <- meta_df$cell[
  meta_df$predictedGroup_Un %in%treg_subtypes  &
    meta_df$Group %in% c("CON", "PP")
]

length(treg_cells)

# 查看每个样本有多少 Treg
table(
  meta_df$Sample[meta_df$cell %in% treg_cells],
  meta_df$Group[meta_df$cell %in% treg_cells]
)

out_treg <- file.path(outDir, "projTreg_only")

unlink(out_treg, recursive = TRUE, force = TRUE)

projTreg <- subsetArchRProject(
  ArchRProj = projHeme2,
  cells = treg_cells,
  outputDirectory = out_treg,
  dropCells = FALSE,
  force = TRUE
)
