library(Seurat)
library(pheatmap)
library(dplyr)
library(RColorBrewer)
TREG_SUBSET <- readRDS("data/TREG_SUBSET.RDS")
DotPlot(TREG_SUBSET,features = genes,group.by = 'treg_sub',cols = c('white','red'))+RotatedAxis()
library(Seurat)
library(pheatmap)
library(dplyr)
FeaturePlot(TREG_SUBSET,features = c('BATF','TNFRSF1B','IL10','SELL'))
genes <- c(
  'BATF','JAK1','STAT3','STAT5','STAT1',
  'IL2RA','TNFRSF9','TNFRSF4','TNFRSF1B','TNFRSF18',
  'SELL','CCR7',
  'FOS','JUN','DUSP1','MAF'
)

# 只保留对象中存在的基因
genes_use <- genes[genes %in% rownames(TREG_SUBSET)]

missing_genes <- setdiff(genes, genes_use)
print(missing_genes)

# 确保分组是你要的亚群
Idents(TREG_SUBSET) <- "treg_sub"

# 计算每个 treg_sub 的平均表达
avg_exp <- AverageExpression(
  TREG_SUBSET,
  features = genes_use,
  group.by = "treg_sub",
  assays = "RNA",
  slot = "data"
)$RNA

# 按基因做 z-score，让热图显示相对高低表达
avg_exp_z <- t(scale(t(avg_exp)))

# 去掉可能产生 NA 的基因
avg_exp_z <- avg_exp_z[complete.cases(avg_exp_z), ]

# 画热图
pheatmap(
  avg_exp_z,
  cluster_rows = TRUE,
  cluster_cols = F,
  
  fontsize_row = 10,
  fontsize_col = 10,
  angle_col = 45,
  border_color = 1,
  main = "Treg subset marker expression"
)

