library(WGCNA)
library(hdWGCNA)
library(tidyverse)
library(cowplot)
library(patchwork)
library(qs)
library(Seurat)
# 基础包
library(Seurat)
library(dplyr)
library(Matrix)

# hdWGCNA
library(hdWGCNA)
library(WGCNA)
####treg####
# 建议：WGCNA 多线程
enableWGCNAThreads()
set.seed(1234)
IBD_T_CELL_remove_DCS <- readRDS("data/IBD_T_CELL_remove_DCS.rds")
seurat_obj <- TREG_SUBSET
# 确认meta列存在
seurat_obj <- SetupForWGCNA(
  seurat_obj,
  gene_select = "fraction", # the gene selection approach
  fraction = 0.05, # fraction of cells that a gene needs to be expressed in order to be included
  wgcna_name = "tutorial" # the name of the hdWGCNA experiment
)
# construct metacells  in each group
seurat_obj <- MetacellsByGroups(
  seurat_obj = seurat_obj,
  group.by = c("T_subtype", "orig.ident"), # specify the columns in seurat_obj@meta.data to group by
  reduction = 'harmony', # select the dimensionality reduction to perform KNN on
  k = 25, # nearest-neighbors parameter
  max_shared = 10, # maximum number of shared cells between two metacells
  ident.group = 'celltype_epi' # set the Idents of the metacell seurat object
)

# normalize metacell expression matrix:
seurat_obj <- NormalizeMetacells(seurat_obj)
seurat_obj <- SetDatExpr(
  seurat_obj,
  group_name = "Treg", # the name of the group of interest in the group.by column
  group.by='T_subtype', # the metadata column containing the cell type info. This same column should have also been used in MetacellsByGroups
  assay = 'RNA', # using RNA assay
  layer = 'data' # using normalized data
)
# Test different soft powers:
seurat_obj <- TestSoftPowers(
  seurat_obj,
  networkType = 'signed' # you can also use "unsigned" or "signed hybrid"
)

# plot the results:
plot_list <- PlotSoftPowers(seurat_obj)

# assemble with patchwork
wrap_plots(plot_list, ncol=2)
power_table <- GetPowerTable(seurat_obj)
head(power_table)
seurat_obj <- ConstructNetwork(
  seurat_obj,
  soft_power = 8, # 自定义了4，如果是自动选择的话可能会是3
  tom_outdir = "TOM",
  tom_name = 'Treg', # name of the topoligical overlap matrix written to disk
  overwrite_tom = TRUE # 允许覆盖已存在的同名文件
)
PlotDendrogram(seurat_obj, main='Treg hdWGCNA Dendrogram')
TOM <- GetTOM(seurat_obj)
seurat_obj <- ModuleEigengenes(
  seurat_obj,
  group.by.vars="orig.ident"
)
# harmonized module eigengenes:
hMEs <- GetMEs(seurat_obj)

# module eigengenes:
MEs <- GetMEs(seurat_obj, harmonized=FALSE)
# compute eigengene-based connectivity (kME):
seurat_obj <- ModuleConnectivity(
  seurat_obj,
  group.by = 'T_subtype', group_name = 'Treg'
)
# rename the modules
seurat_obj <- ResetModuleNames(
  seurat_obj,
  new_name = "Treg"
)
# plot genes ranked by kME for each module
p <- PlotKMEs(seurat_obj, ncol=5)

p
# get the module assignment table:
modules <- GetModules(seurat_obj) %>% subset(module != 'grey')

# show the first 6 columns:
head(modules[,1:6])
# get hub genes
hub_df <- GetHubGenes(seurat_obj, n_hubs = 10)

head(hub_df)
# compute gene scoring for the top 25 hub genes by kME for each module
# with UCell method
library(UCell)
seurat_obj <- ModuleExprScore(
  seurat_obj,
  n_genes = 25,
  method='UCell'
)
saveRDS(seurat_obj, file='hdWGCNA_object.rds')
ModuleCorrelogram(seurat_obj)
###RESTART####
hdWGCNA_object_TREG <- readRDS("data/hdWGCNA_object_TREG.rds")
seurat_obj <- hdWGCNA_object_TREG

# 每个模块制作hMEs的特征图
plot_list <- ModuleFeaturePlot(
  seurat_obj,
  features='hMEs', # plot the hMEs
  order=TRUE # order so the points with highest hMEs are on top
)

# stitch together with patchwork
wrap_plots(plot_list, ncol=4)

# 每个模块制作hub scores的特征图
plot_list <- ModuleFeaturePlot(
  seurat_obj,
  features='scores', # plot the hub gene scores
  order='shuffle', # order so cells are shuffled
  ucell = TRUE # depending on Seurat vs UCell for gene scoring
)
# stitch together with patchwork
wrap_plots(plot_list, ncol=4)

# 每个模块在不同样本中的情况
seurat_obj$cluster <- do.call(rbind, strsplit(as.character(seurat_obj$orig.ident), ' '))[,1]

ModuleRadarPlot(
  seurat_obj,
  group.by = 'cluster',
  barcodes = seurat_obj@meta.data %>% 
    subset(celltype == 'Treg') %>% 
    rownames(),
  axis.label.size=4,
  grid.label.size=4
)

# 查看模块相关图
ModuleCorrelogram(seurat_obj)

# get hMEs from seurat object
MEs <- GetMEs(seurat_obj, harmonized=TRUE)
modules <- GetModules(seurat_obj)
mods <- levels(modules$module); mods <- mods[mods != 'grey']

# add hMEs to Seurat meta-data:
seurat_obj@meta.data <- cbind(seurat_obj@meta.data, MEs)

# plot with Seurat's DotPlot function
treg <- subset(seurat_obj,T_subtype=='Treg')
treg$cells <- rownames(treg@meta.data)
treg@meta.data <- treg@meta.data %>% left_join(y = TREG_SUBSET@meta.data,by = 'cell_name')
rownames(treg@meta.data) <- treg$cells

p <- DotPlot(treg, features=mods, group.by = 'treg_sub')

# flip the x/y axes, rotate the axis labels, and change color scheme:
p <- p +
  RotatedAxis() +
  scale_color_gradient2(high='red', mid='grey95', low='blue')

# plot output
p

# 使用ModuleNetworkPlot可视化每个模块前50(数值可自定)的hub gene
ModuleNetworkPlot(
  seurat_obj, 
  outdir='ModuleNetworks', # new folder name
  n_inner = 20, # number of genes in inner ring
  n_outer = 30, # number of genes in outer ring
  n_conns = Inf, # show all of the connections
  plot_size=c(10,10), # larger plotting area
  vertex.label.cex=1 # font size
)

options(future.globals.maxSize = 5 * 1024^3)  # 5GB

# hubgene network(基因数可自定)
HubGeneNetworkPlot(
  seurat_obj,
  n_hubs = 2, 
  n_other=2,
  edge_prop = 0.75,
  mods = 'all'
)

# 可以选择模块数
g <- HubGeneNetworkPlot(seurat_obj,  return_graph=TRUE)
# get the list of modules:
modules <- GetModules(seurat_obj)
mods <- levels(modules$module); mods <- mods[mods != 'grey']
# hubgene network
HubGeneNetworkPlot(
  seurat_obj,
  n_hubs = 2, 
  n_other= 2,
  edge_prop = 0.75,
  mods = mods[1:5] # only select 5 modules
)

seurat_obj <- RunModuleUMAP(
  seurat_obj,
  n_hubs = 10, # number of hub genes to include for the UMAP embedding
  n_neighbors=15, # neighbors parameter for UMAP
  min_dist=0.1 # min distance between points in UMAP space
)

# get the hub gene UMAP table from the seurat object
umap_df <- GetModuleUMAP(seurat_obj)

# plot with ggplot
ggplot(umap_df, aes(x=UMAP1, y=UMAP2)) +
  geom_point(
    color=umap_df$color, # color each point by WGCNA module
    size=umap_df$kME*2 # size of each point based on intramodular connectivity
  ) +
  umap_theme()


ModuleUMAPPlot(
  seurat_obj,
  edge.alpha=0.25,
  sample_edges=TRUE,
  edge_prop=0.1, # proportion of edges to sample (20% here)
  label_hubs=2 ,# how many hub genes to plot per module?
  keep_grey_edges=FALSE
)


library(Seurat)
library(tidyverse)
library(cowplot)
library(patchwork)
library(WGCNA)
library(hdWGCNA)
library(enrichR)
library(GeneOverlap)
library(qs)

#dir.create("14-hdWGCNA")
#setwd("14-hdWGCNA")

seurat_obj <- qread("hdWGCNA_object.qs")

# 定义enrichr databases
dbs <- c('GO_Biological_Process_2023',
         'GO_Cellular_Component_2023',
         'GO_Molecular_Function_2023')

# 富集分析
seurat_obj <- RunEnrichr(
  seurat_obj,
  dbs=dbs,
  max_genes = 100 # use max_genes = Inf to choose all genes
)

# 检索输出表
enrich_df <- GetEnrichrTable(seurat_obj)

# 查看结果
head(enrich_df)

# make GO term plots:
EnrichrBarPlot(
  seurat_obj,
  outdir = "enrichr_plots", # name of output directory
  n_terms = 10, # number of enriched terms to show (sometimes more are shown if there are ties)
  plot_size = c(5,7), # width, height of the output .pdfs
  logscale=TRUE # do you want to show the enrichment as a log scale?
)

# enrichr dotplot
EnrichrDotPlot(
  seurat_obj,
  mods = "all", # use all modules (default)
  database = "GO_Biological_Process_2023", # this must match one of the dbs used previously
  n_terms=2, # number of terms per module
  term_size=8, # font size for the terms
  p_adj = FALSE # show the p-val or adjusted p-val?
)  + scale_color_stepsn(colors=rev(viridis::magma(256)))














####判断其后炎症的相关性####
seurat_obj$inflammation_score <- seurat_obj$characteristics_ch1.14

# 检查是否有对不上的
sum(is.na(seurat_obj$inflammation_score))
# 如 NA 很多，说明两个对象的 cell barcode 不一致，需要先统一/映射

# 确保是数值型

## 1) 计算模块-表型相关（按细胞亚群分别算；你这里用 T_subtype，能单独得到 Treg 的相关矩阵）
seurat_obj_meta <- seurat_obj@meta.data%>%
  mutate(
    inflammation_score = as.character(inflammation_score),
    # 如果里面可能有 "score: 3.2" 这种，建议用 parse_number；纯数字用 as.numeric 也行
    inflammation_score = readr::parse_number(inflammation_score)
    # inflammation_score = as.numeric(inflammation_score)  # 纯数字字符串时也可
  ) 
seurat_obj@meta.data <- seurat_obj_meta
## 0) 保证 traits 是 numeric（NA 保留没关系）
x <- as.character(seurat_obj$inflammation_score)
x[x %in% c("NA","NaN","Inf","-Inf","", "NULL")] <- NA
x <- gsub(",", "", x)
seurat_obj$inflammation_score <- as.numeric(x)
seurat_obj$inflammation_score[!is.finite(seurat_obj$inflammation_score)] <- NA

## 1) group.by 必须是 factor（建议显式转一下）
seurat_obj$T_subtype <- as.factor(seurat_obj$T_subtype)

## 2) 做一个“是否有评分”的筛选列
seurat_obj$has_inflam_score <- factor(
  ifelse(is.na(seurat_obj$inflammation_score), "no", "yes"),
  levels = c("no", "yes")
)
# 确保是 numeric
seurat_obj$inflammation_score <- as.numeric(as.character(seurat_obj$inflammation_score))
seurat_obj$inflammation_score[!is.finite(seurat_obj$inflammation_score)] <- NA

# 用函数内部 subset（不要先 subset()）
seurat_obj <- ModuleTraitCorrelation(
  seurat_obj,
  traits = "inflammation_score",
  group.by = "T_subtype",
  subset_by = "has_inflam_score",
  subset_groups = "yes"
)










# 取 harmonized MEs（metacell层面通常用这个）
MEs <- GetMEs(seurat_obj, harmonized = TRUE)

# MEs 的行名应该是细胞/ metacell 的名字；把它拼到 meta.data
md <- seurat_obj@meta.data

# 对齐（取交集，避免不匹配）
common <- intersect(rownames(md), rownames(MEs))
md2 <- md[common, , drop=FALSE]
MEs2 <- MEs[common, , drop=FALSE]

# 只保留有炎症评分的
md2$inflammation_score <- as.numeric(as.character(md2$inflammation_score))
keep <- is.finite(md2$inflammation_score)
md2 <- md2[keep, , drop=FALSE]
MEs2 <- MEs2[keep, , drop=FALSE]


library(dplyr)
library(purrr)
library(tibble)

calc_group <- function(idx){
  xmat <- MEs2[idx, , drop=FALSE]
  y <- md2$inflammation_score[idx]
  map_dfr(colnames(xmat), \(m){
    x <- xmat[, m]
    ok <- is.finite(x) & is.finite(y)
    n <- sum(ok)
    if(n < 3) return(tibble(module=m, n=n, cor=NA_real_, p=NA_real_))
    ct <- suppressWarnings(cor.test(x[ok], y[ok], method="spearman", exact=FALSE))
    tibble(module=m, n=n, cor=unname(ct$estimate), p=ct$p.value)
  })
}

res <- md2 %>%
  mutate(.cell = rownames(md2)) %>%
  group_by(T_subtype) %>%
  group_modify(\(df, key){
    idx <- match(df$.cell, rownames(MEs2))
    out <- calc_group(idx)
    out
  }) %>%
  ungroup() %>%
  group_by(T_subtype) %>%
  mutate(fdr = p.adjust(p, method="fdr")) %>%
  ungroup()


library(ggplot2)

res2 <- res %>%
  filter(module != "grey") %>%
  mutate(star = case_when(
    fdr < 0.001 ~ "***",
    fdr < 0.01 ~ "**",
    fdr < 0.05 ~ "*",
    TRUE ~ ""
  ))

p <- ggplot(res2, aes(x=T_subtype, y=module, fill=cor)) +
  geom_tile() +
  geom_text(aes(label=star), size=3) +
  theme_bw() +
  theme(axis.text.x = element_text(angle=45, hjust=1)) +
  labs(x="T_subtype", y="Module (hME)", fill="Spearman r",
       title="Module vs inflammation_score (only cells with score)")
setwd("data\\UC")  
p
ggsave("Module_inflammation_tile.pdf", p, width=10, height=8)
library(dplyr)

# 1) 模块基因表（包含每个基因的 module、kME 等）
mods <- GetModules(seurat_obj)
write.csv(mods, "treg_top_hubs.csv", row.names = FALSE)
# 看看模块名到底是 "treg5" 还是 "Treg5"
levels(mods$module)

# 2) 取 treg5 的基因，并按 kME 从大到小排
treg3_genes <- mods %>%
  filter(tolower(module) == "treg1") %>%   # 不区分大小写匹配
  arrange(desc(kME_Treg1)) %>%
  mutate(rank = row_number())

# 3) hub 基因：取前 N 个
n_hub <- 20
treg1_hubs <- treg3_genes %>% slice_head(n = n_hub)

treg1_hubs %>% select(gene_name, module, module, rank) %>% head(20)

# 导出
write.csv(treg5_hubs, "treg5_top_hubs.csv", row.names = FALSE)


####判断其和疾病的相关性####


seurat_obj1 <- ModuleTraitCorrelation(
  seurat_obj,
  traits = "disease",
  group.by = "T_subtype"
)
# 转 factor
seurat_obj$disease <- as.factor(seurat_obj$disease)

# 再跑
seurat_obj1 <- ModuleTraitCorrelation(
  seurat_obj,
  traits = "disease",
  group.by = "T_subtype"
)
MEs <- GetMEs(seurat_obj1, harmonized = TRUE)  # 建议用 harmonized
md  <- seurat_obj@meta.data

common <- intersect(rownames(md), rownames(MEs))
md2 <- md[common, , drop=FALSE]
MEs2 <- MEs[common, , drop=FALSE]

df <- cbind(md2, MEs2)
mod_cols <- colnames(MEs2)
head(mod_cols, 10)
library(ggplot2)

# 找到对应列名（按你实际列名调整）
treg1_col <- grep("Treg1", mod_cols, ignore.case = TRUE, value = TRUE)[1]

p_vln <- ggplot(df, aes(x = disease, y = .data[[treg1_col]])) +
  geom_violin(trim=FALSE) +
  geom_boxplot(width=0.15, outlier.size=0.2) +
  facet_wrap(~T_subtype, scales="free_y") +
  theme_bw() +
  theme(axis.text.x = element_text(angle=45, hjust=1)) +
  labs(x="Disease", y=treg1_col, title="Module activity by disease (facet by T_subtype)")

p_vln
ggsave("disease_Treg1_module_violin.pdf", p_vln, width=12, height=7)
library(dplyr)
library(tidyr)
library(purrr)
library(tibble)

# 确保是二分类（多分类可以先选两组，或做ANOVA/Kruskal）
levels(df$disease)

# 你想比较的两组（按你的实际标签改）
g1 <- levels(df$disease)[1]
g2 <- levels(df$disease)[2]

calc_effect <- function(d, m){
  x <- d[[m]]
  grp <- d$disease
  ok <- is.finite(x) & !is.na(grp)
  d <- d[ok, , drop=FALSE]
  if(length(unique(d$disease)) < 2) return(tibble(module=m, n=0, effect=NA_real_, p=NA_real_))
  
  # 用中位数差作为效应量（g2 - g1）
  med1 <- median(d[[m]][d$disease == g1])
  med2 <- median(d[[m]][d$disease == g2])
  eff  <- med2 - med1
  
  p <- suppressWarnings(wilcox.test(d[[m]] ~ d$disease)$p.value)
  tibble(module=m, n=nrow(d), effect=eff, p=p)
}

res <- df %>%
  group_by(T_subtype) %>%
  group_modify(~map_dfr(mod_cols, \(m) calc_effect(.x, m))) %>%
  ungroup() %>%
  group_by(T_subtype) %>%
  mutate(fdr = p.adjust(p, method="fdr"),
         star = case_when(fdr < 0.001 ~ "***",
                          fdr < 0.01 ~ "**",
                          fdr < 0.05 ~ "*",
                          TRUE ~ "")) %>%
  ungroup()

p_tile <- ggplot(res, aes(x=T_subtype, y=module, fill=effect)) +
  geom_tile() +
  geom_text(aes(label=star), size=3) +
  theme_bw() +
  theme(axis.text.x = element_text(angle=45, hjust=1)) +
  labs(x="T_subtype", y="Module (hME)", fill=paste0("Median(", g2, ") - Median(", g1, ")"),
       title="Disease effect on module activity (per T_subtype)")

p_tile
ggsave("disease_module_effect_tile.pdf", p_tile, width=12, height=0.25*length(mod_cols)+3)
library(WGCNA)
library(hdWGCNA)
library(tidyverse)
library(cowplot)
library(patchwork)
library(qs)
library(Seurat)
# 基础包
library(Seurat)
library(dplyr)
library(Matrix)

# hdWGCNA
library(hdWGCNA)
library(WGCNA)
####TH17####
# 建议：WGCNA 多线程
enableWGCNAThreads()
set.seed(1234)
IBD_T_CELL_remove_DCS <- readRDS("data/IBD_T_CELL_remove_DCS.rds")
seurat_obj <- sce.all.filt
# 确认meta列存在
seurat_obj <- SetupForWGCNA(
  seurat_obj,
  gene_select = "fraction", # the gene selection approach
  fraction = 0.05, # fraction of cells that a gene needs to be expressed in order to be included
  wgcna_name = "tutorial" # the name of the hdWGCNA experiment
)
# construct metacells  in each group
seurat_obj <- MetacellsByGroups(
  seurat_obj = seurat_obj,
  group.by = c("T_subtype", "orig.ident"), # specify the columns in seurat_obj@meta.data to group by
  reduction = 'harmony', # select the dimensionality reduction to perform KNN on
  k = 25, # nearest-neighbors parameter
  max_shared = 10, # maximum number of shared cells between two metacells
  ident.group = 'T_subtype' # set the Idents of the metacell seurat object
)

# normalize metacell expression matrix:
seurat_obj <- NormalizeMetacells(seurat_obj)
seurat_obj <- SetDatExpr(
  seurat_obj,
  group_name = "TH17", # the name of the group of interest in the group.by column
  group.by='T_subtype', # the metadata column containing the cell type info. This same column should have also been used in MetacellsByGroups
  assay = 'RNA', # using RNA assay
  layer = 'data' # using normalized data
)
# Test different soft powers:
seurat_obj <- TestSoftPowers(
  seurat_obj,
  networkType = 'signed' # you can also use "unsigned" or "signed hybrid"
)

# plot the results:
plot_list <- PlotSoftPowers(seurat_obj)

# assemble with patchwork
wrap_plots(plot_list, ncol=2)
power_table <- GetPowerTable(seurat_obj)
head(power_table)
seurat_obj <- ConstructNetwork(
  seurat_obj,
  soft_power = 9, # 自定义了4，如果是自动选择的话可能会是3
  tom_outdir = "TOM",
  tom_name = 'Treg', # name of the topoligical overlap matrix written to disk
  overwrite_tom = TRUE # 允许覆盖已存在的同名文件
)
PlotDendrogram(seurat_obj, main='Treg hdWGCNA Dendrogram')
TOM <- GetTOM(seurat_obj)
seurat_obj <- ModuleEigengenes(
  seurat_obj,
  group.by.vars="orig.ident"
)
# harmonized module eigengenes:
hMEs <- GetMEs(seurat_obj)

# module eigengenes:
MEs <- GetMEs(seurat_obj, harmonized=FALSE)
# compute eigengene-based connectivity (kME):
seurat_obj <- ModuleConnectivity(
  seurat_obj,
  group.by = 'T_subtype', group_name = 'TH17'
)
# rename the modules
seurat_obj <- ResetModuleNames(
  seurat_obj,
  new_name = "TH17"
)
# plot genes ranked by kME for each module
p <- PlotKMEs(seurat_obj, ncol=5)

p
# get the module assignment table:
modules <- GetModules(seurat_obj) %>% subset(module != 'grey')

# show the first 6 columns:
head(modules[,1:6])
# get hub genes
hub_df <- GetHubGenes(seurat_obj, n_hubs = 10)

head(hub_df)
# compute gene scoring for the top 25 hub genes by kME for each module
# with UCell method
library(UCell)
seurat_obj <- ModuleExprScore(
  seurat_obj,
  n_genes = 25,
  method='UCell'
)
saveRDS(seurat_obj, file='hdWGCNA_object_TH17.rds')

seurat_obj <- hdWGCNA_object

####判断其后炎症的相关性####
seurat_obj$inflammation_score <- seurat_obj$characteristics_ch1.14

# 检查是否有对不上的
sum(is.na(seurat_obj$inflammation_score))
# 如 NA 很多，说明两个对象的 cell barcode 不一致，需要先统一/映射

# 确保是数值型

## 1) 计算模块-表型相关（按细胞亚群分别算；你这里用 T_subtype，能单独得到 Treg 的相关矩阵）
seurat_obj_meta <- seurat_obj@meta.data%>%
  mutate(
    inflammation_score = as.character(inflammation_score),
    # 如果里面可能有 "score: 3.2" 这种，建议用 parse_number；纯数字用 as.numeric 也行
    inflammation_score = readr::parse_number(inflammation_score)
    # inflammation_score = as.numeric(inflammation_score)  # 纯数字字符串时也可
  ) 
seurat_obj@meta.data <- seurat_obj_meta
## 0) 保证 traits 是 numeric（NA 保留没关系）
x <- as.character(seurat_obj$inflammation_score)
x[x %in% c("NA","NaN","Inf","-Inf","", "NULL")] <- NA
x <- gsub(",", "", x)
seurat_obj$inflammation_score <- as.numeric(x)
seurat_obj$inflammation_score[!is.finite(seurat_obj$inflammation_score)] <- NA

## 1) group.by 必须是 factor（建议显式转一下）
seurat_obj$T_subtype <- as.factor(seurat_obj$T_subtype)

## 2) 做一个“是否有评分”的筛选列
seurat_obj$has_inflam_score <- factor(
  ifelse(is.na(seurat_obj$inflammation_score), "no", "yes"),
  levels = c("no", "yes")
)
# 确保是 numeric
seurat_obj$inflammation_score <- as.numeric(as.character(seurat_obj$inflammation_score))
seurat_obj$inflammation_score[!is.finite(seurat_obj$inflammation_score)] <- NA

# 用函数内部 subset（不要先 subset()）
seurat_obj <- ModuleTraitCorrelation(
  seurat_obj,
  traits = "inflammation_score",
  group.by = "T_subtype",
  subset_by = "has_inflam_score",
  subset_groups = "yes"
)










# 取 harmonized MEs（metacell层面通常用这个）
MEs <- GetMEs(seurat_obj, harmonized = TRUE)

# MEs 的行名应该是细胞/ metacell 的名字；把它拼到 meta.data
md <- seurat_obj@meta.data

# 对齐（取交集，避免不匹配）
common <- intersect(rownames(md), rownames(MEs))
md2 <- md[common, , drop=FALSE]
MEs2 <- MEs[common, , drop=FALSE]

# 只保留有炎症评分的
md2$inflammation_score <- as.numeric(as.character(md2$inflammation_score))
keep <- is.finite(md2$inflammation_score)
md2 <- md2[keep, , drop=FALSE]
MEs2 <- MEs2[keep, , drop=FALSE]


library(dplyr)
library(purrr)
library(tibble)

calc_group <- function(idx){
  xmat <- MEs2[idx, , drop=FALSE]
  y <- md2$inflammation_score[idx]
  map_dfr(colnames(xmat), \(m){
    x <- xmat[, m]
    ok <- is.finite(x) & is.finite(y)
    n <- sum(ok)
    if(n < 3) return(tibble(module=m, n=n, cor=NA_real_, p=NA_real_))
    ct <- suppressWarnings(cor.test(x[ok], y[ok], method="spearman", exact=FALSE))
    tibble(module=m, n=n, cor=unname(ct$estimate), p=ct$p.value)
  })
}

res <- md2 %>%
  mutate(.cell = rownames(md2)) %>%
  group_by(T_subtype) %>%
  group_modify(\(df, key){
    idx <- match(df$.cell, rownames(MEs2))
    out <- calc_group(idx)
    out
  }) %>%
  ungroup() %>%
  group_by(T_subtype) %>%
  mutate(fdr = p.adjust(p, method="fdr")) %>%
  ungroup()


library(ggplot2)

res2 <- res %>%
  filter(module != "grey") %>%
  mutate(star = case_when(
    fdr < 0.001 ~ "***",
    fdr < 0.01 ~ "**",
    fdr < 0.05 ~ "*",
    TRUE ~ ""
  ))

p <- ggplot(res2, aes(x=T_subtype, y=module, fill=cor)) +
  geom_tile() +
  geom_text(aes(label=star), size=3) +
  theme_bw() +
  theme(axis.text.x = element_text(angle=45, hjust=1)) +
  labs(x="T_subtype", y="Module (hME)", fill="Spearman r",
       title="Module vs inflammation_score (only cells with score)")
setwd("data\\UC")  
p
ggsave("Module_inflammation_tile.pdf", p, width=10, height=8)
library(dplyr)

# 1) 模块基因表（包含每个基因的 module、kME 等）
mods <- GetModules(seurat_obj)
write.csv(mods, "TH17_genes_top_hubs.csv", row.names = FALSE)
# 看看模块名到底是 "treg5" 还是 "Treg5"
levels(mods$module)

# 2) 取 treg5 的基因，并按 kME 从大到小排
TH175_genes <- mods %>%
  filter(module== "TH175") %>%   # 不区分大小写匹配
  arrange(desc(kME_TH175)) %>%
  mutate(rank = row_number())

# 3) hub 基因：取前 N 个
n_hub <- 20
TH175_genes <- TH175_genes %>% slice_head(n = n_hub)

TH175_genes %>% select(gene_name, module, module, rank) %>% head(20)

# 导出
write.csv(TH175_genes, "TH175_genes_top_hubs.csv", row.names = FALSE)
TH173_genes <- mods %>%
  filter(module== "TH173") %>%   # 不区分大小写匹配
  arrange(desc(kME_TH173)) %>%
  mutate(rank = row_number())

# 3) hub 基因：取前 N 个
n_hub <- 20
TH173_genes <- TH173_genes %>% slice_head(n = n_hub)

TH173_genes %>% select(gene_name, module, module, rank) %>% head(20)
write.csv(TH173_genes, "TH173_genes_top_hubs.csv", row.names = FALSE)





####判断其和疾病的相关性####


seurat_obj1 <- ModuleTraitCorrelation(
  seurat_obj,
  traits = "disease",
  group.by = "T_subtype"
)
# 转 factor
seurat_obj$disease <- as.factor(seurat_obj$disease)

# 再跑
seurat_obj1 <- ModuleTraitCorrelation(
  seurat_obj,
  traits = "disease",
  group.by = "T_subtype"
)
MEs <- GetMEs(seurat_obj1, harmonized = TRUE)  # 建议用 harmonized
md  <- seurat_obj@meta.data

common <- intersect(rownames(md), rownames(MEs))
md2 <- md[common, , drop=FALSE]
MEs2 <- MEs[common, , drop=FALSE]

df <- cbind(md2, MEs2)
mod_cols <- colnames(MEs2)
head(mod_cols, 10)
library(ggplot2)

# 找到对应列名（按你实际列名调整）
treg1_col <- grep("Treg1", mod_cols, ignore.case = TRUE, value = TRUE)[1]

p_vln <- ggplot(df, aes(x = disease, y = .data[[treg1_col]])) +
  geom_violin(trim=FALSE) +
  geom_boxplot(width=0.15, outlier.size=0.2) +
  facet_wrap(~T_subtype, scales="free_y") +
  theme_bw() +
  theme(axis.text.x = element_text(angle=45, hjust=1)) +
  labs(x="Disease", y=treg1_col, title="Module activity by disease (facet by T_subtype)")

p_vln
ggsave("disease_Treg1_module_violin.pdf", p_vln, width=12, height=7)
library(dplyr)
library(tidyr)
library(purrr)
library(tibble)

# 确保是二分类（多分类可以先选两组，或做ANOVA/Kruskal）
levels(df$disease)

# 你想比较的两组（按你的实际标签改）
g1 <- levels(df$disease)[1]
g2 <- levels(df$disease)[2]

calc_effect <- function(d, m){
  x <- d[[m]]
  grp <- d$disease
  ok <- is.finite(x) & !is.na(grp)
  d <- d[ok, , drop=FALSE]
  if(length(unique(d$disease)) < 2) return(tibble(module=m, n=0, effect=NA_real_, p=NA_real_))
  
  # 用中位数差作为效应量（g2 - g1）
  med1 <- median(d[[m]][d$disease == g1])
  med2 <- median(d[[m]][d$disease == g2])
  eff  <- med2 - med1
  
  p <- suppressWarnings(wilcox.test(d[[m]] ~ d$disease)$p.value)
  tibble(module=m, n=nrow(d), effect=eff, p=p)
}

res <- df %>%
  group_by(T_subtype) %>%
  group_modify(~map_dfr(mod_cols, \(m) calc_effect(.x, m))) %>%
  ungroup() %>%
  group_by(T_subtype) %>%
  mutate(fdr = p.adjust(p, method="fdr"),
         star = case_when(fdr < 0.001 ~ "***",
                          fdr < 0.01 ~ "**",
                          fdr < 0.05 ~ "*",
                          TRUE ~ "")) %>%
  ungroup()

p_tile <- ggplot(res, aes(x=T_subtype, y=module, fill=effect)) +
  geom_tile() +
  geom_text(aes(label=star), size=3) +
  theme_bw() +
  theme(axis.text.x = element_text(angle=45, hjust=1)) +
  labs(x="T_subtype", y="Module (hME)", fill=paste0("Median(", g2, ") - Median(", g1, ")"),
       title="Disease effect on module activity (per T_subtype)")

p_tile
ggsave("disease_module_effect_tile.pdf", p_tile, width=12, height=0.25*length(mod_cols)+3)



###取交集

Treg1_genes <- treg_top_hubs %>%
  filter(module== "Treg1") %>%   # 不区分大小写匹配
  arrange(desc(kME_Treg1)) %>%
  mutate(rank = row_number())

# 3) hub 基因：取前 N 个
n_hub <- 30
Treg1_genes <- Treg1_genes %>% slice_head(n = n_hub)




Treg2_genes <- treg_top_hubs %>%
  filter(module== "Treg2") %>%   # 不区分大小写匹配
  arrange(desc(kME_Treg2)) %>%
  mutate(rank = row_number())

# 3) hub 基因：取前 N 个
n_hub <- 30
Treg2_genes <- Treg2_genes %>% slice_head(n = n_hub)


Treg5_genes <- treg_top_hubs %>%
  filter(module== "Treg5") %>%   # 不区分大小写匹配
  arrange(desc(kME_Treg5)) %>%
  mutate(rank = row_number())

# 3) hub 基因：取前 N 个
n_hub <- 30
Treg5_genes <- Treg5_genes %>% slice_head(n = n_hub)

Treg6_genes <- treg_top_hubs %>%
  filter(module== "Treg6") %>%   # 不区分大小写匹配
  arrange(desc(kME_Treg6)) %>%
  mutate(rank = row_number())

# 3) hub 基因：取前 N 个
n_hub <- 30
Treg6_genes <- Treg6_genes %>% slice_head(n = n_hub)





TH172_genes <- TH17_genes_top_hubs %>%
  filter(module== "TH172") %>%   # 不区分大小写匹配
  arrange(desc(kME_TH172)) %>%
  mutate(rank = row_number())

# 3) hub 基因：取前 N 个
n_hub <- 30
TH172_genes <- TH172_genes %>% slice_head(n = n_hub)


TH175_genes <- TH17_genes_top_hubs %>%
  filter(module== "TH175") %>%   # 不区分大小写匹配
  arrange(desc(kME_TH175)) %>%
  mutate(rank = row_number())

# 3) hub 基因：取前 N 个
n_hub <- 30
TH175_genes <- TH175_genes %>% slice_head(n = n_hub)

intersect(Treg1_genes$gene_name,TH172_genes$gene_name)
intersect(Treg1_genes$gene_name,TH175_genes$gene_name)
co_factors <- c(intersect(Treg1_genes$gene_name,TH172_genes$gene_name),intersect(Treg1_genes$gene_name,TH175_genes$gene_name))
IBD_related_gene <- read_excel("data/IBD相关基因集-分析(2).xlsx",sheet = "数据清洗")
head(IBD_related_gene)
gene_vector2 <- unlist(IBD_related_gene, use.names = FALSE)
intersect(co_factors,gene_vector2)
IBD_T_CELL_remove_dcs2$CELL_DISE <- paste0(IBD_T_CELL_remove_dcs2$T_subtype,'_',IBD_T_CELL_remove_dcs2$disease)
DotPlot(IBD_T_CELL_remove_dcs2,features = co_factors,group.by = 'CELL_DISE')


library(WGCNA)
library(hdWGCNA)
library(tidyverse)
library(cowplot)
library(patchwork)
library(qs)
library(Seurat)
# 基础包
library(Seurat)
library(dplyr)
library(Matrix)

# hdWGCNA
library(hdWGCNA)
library(WGCNA)
####trm####
# 建议：WGCNA 多线程
enableWGCNAThreads()
set.seed(1234)
IBD_T_CELL_remove_dcs2 <- readRDS("data/IBD_T_CELL_remove_dcs2.rds")
seurat_obj <- IBD_T_CELL_remove_dcs2
# 确认meta列存在
seurat_obj <- SetupForWGCNA(
  seurat_obj,
  gene_select = "fraction", # the gene selection approach
  fraction = 0.05, # fraction of cells that a gene needs to be expressed in order to be included
  wgcna_name = "tutorial" # the name of the hdWGCNA experiment
)
# construct metacells  in each group
seurat_obj <- MetacellsByGroups(
  seurat_obj = seurat_obj,
  group.by = c("T_subtype", "orig.ident"), # specify the columns in seurat_obj@meta.data to group by
  reduction = 'harmony', # select the dimensionality reduction to perform KNN on
  k = 25, # nearest-neighbors parameter
  max_shared = 10, # maximum number of shared cells between two metacells
  ident.group = 'T_subtype' # set the Idents of the metacell seurat object
)

# normalize metacell expression matrix:
seurat_obj <- NormalizeMetacells(seurat_obj)
seurat_obj <- SetDatExpr(
  seurat_obj,
  group_name = "Trm_CD8", # the name of the group of interest in the group.by column
  group.by='T_subtype', # the metadata column containing the cell type info. This same column should have also been used in MetacellsByGroups
  assay = 'RNA', # using RNA assay
  layer = 'data' # using normalized data
)
library(future)
plan(sequential)
library(WGCNA)
disableWGCNAThreads()
# 或者
allowWGCNAThreads(nThreads = 1)

# Test different soft powers:
seurat_obj <- TestSoftPowers(
  seurat_obj,
  networkType = 'signed' # you can also use "unsigned" or "signed hybrid"
)

# plot the results:
plot_list <- PlotSoftPowers(seurat_obj)

# assemble with patchwork
wrap_plots(plot_list, ncol=2)
power_table <- GetPowerTable(seurat_obj)
head(power_table)
seurat_obj <- ConstructNetwork(
  seurat_obj,
  soft_power = 7, # 自定义了4，如果是自动选择的话可能会是3
  tom_outdir = "TOM",
  tom_name = 'Trm_CD8', # name of the topoligical overlap matrix written to disk
  overwrite_tom = TRUE # 允许覆盖已存在的同名文件
)
PlotDendrogram(seurat_obj, main='Trm_CD8 hdWGCNA Dendrogram')
TOM <- GetTOM(seurat_obj)
seurat_obj <- ModuleEigengenes(
  seurat_obj,
  group.by.vars="orig.ident"
)
# harmonized module eigengenes:
hMEs <- GetMEs(seurat_obj)

# module eigengenes:
MEs <- GetMEs(seurat_obj, harmonized=FALSE)
# compute eigengene-based connectivity (kME):
seurat_obj <- ModuleConnectivity(
  seurat_obj,
  group.by = 'T_subtype', group_name = 'Trm_CD8'
)
# rename the modules
seurat_obj <- ResetModuleNames(
  seurat_obj,
  new_name = "Trm_CD8_"
)
# plot genes ranked by kME for each module
p <- PlotKMEs(seurat_obj, ncol=5)

p
# get the module assignment table:
modules <- GetModules(seurat_obj) %>% subset(module != 'grey')

# show the first 6 columns:
head(modules[,1:6])
# get hub genes
hub_df <- GetHubGenes(seurat_obj, n_hubs = 10)

head(hub_df)
# compute gene scoring for the top 25 hub genes by kME for each module
# with UCell method
library(UCell)
seurat_obj <- ModuleExprScore(
  seurat_obj,
  n_genes = 25,
  method='UCell'
)
saveRDS(seurat_obj, file='hdWGCNA_object_trmcd8.rds')

#seurat_obj <- hdWGCNA_object_TREG

####判断其后炎症的相关性####
seurat_obj$inflammation_score <- seurat_obj$characteristics_ch1.14

# 检查是否有对不上的
sum(is.na(seurat_obj$inflammation_score))
# 如 NA 很多，说明两个对象的 cell barcode 不一致，需要先统一/映射

# 确保是数值型

## 1) 计算模块-表型相关（按细胞亚群分别算；你这里用 T_subtype，能单独得到 Treg 的相关矩阵）
seurat_obj_meta <- seurat_obj@meta.data%>%
  mutate(
    inflammation_score = as.character(inflammation_score),
    # 如果里面可能有 "score: 3.2" 这种，建议用 parse_number；纯数字用 as.numeric 也行
    inflammation_score = readr::parse_number(inflammation_score)
    # inflammation_score = as.numeric(inflammation_score)  # 纯数字字符串时也可
  ) 
seurat_obj@meta.data <- seurat_obj_meta
## 0) 保证 traits 是 numeric（NA 保留没关系）
x <- as.character(seurat_obj$inflammation_score)
x[x %in% c("NA","NaN","Inf","-Inf","", "NULL")] <- NA
x <- gsub(",", "", x)
seurat_obj$inflammation_score <- as.numeric(x)
seurat_obj$inflammation_score[!is.finite(seurat_obj$inflammation_score)] <- NA

## 1) group.by 必须是 factor（建议显式转一下）
seurat_obj$T_subtype <- as.factor(seurat_obj$T_subtype)

## 2) 做一个“是否有评分”的筛选列
seurat_obj$has_inflam_score <- factor(
  ifelse(is.na(seurat_obj$inflammation_score), "no", "yes"),
  levels = c("no", "yes")
)
# 确保是 numeric
seurat_obj$inflammation_score <- as.numeric(as.character(seurat_obj$inflammation_score))
seurat_obj$inflammation_score[!is.finite(seurat_obj$inflammation_score)] <- NA

# 用函数内部 subset（不要先 subset()）
seurat_obj <- ModuleTraitCorrelation(
  seurat_obj,
  traits = "inflammation_score",
  group.by = "T_subtype",
  subset_by = "has_inflam_score",
  subset_groups = "yes"
)










# 取 harmonized MEs（metacell层面通常用这个）
MEs <- GetMEs(seurat_obj, harmonized = TRUE)

# MEs 的行名应该是细胞/ metacell 的名字；把它拼到 meta.data
md <- seurat_obj@meta.data

# 对齐（取交集，避免不匹配）
common <- intersect(rownames(md), rownames(MEs))
md2 <- md[common, , drop=FALSE]
MEs2 <- MEs[common, , drop=FALSE]

# 只保留有炎症评分的
md2$inflammation_score <- as.numeric(as.character(md2$inflammation_score))
keep <- is.finite(md2$inflammation_score)
md2 <- md2[keep, , drop=FALSE]
MEs2 <- MEs2[keep, , drop=FALSE]


library(dplyr)
library(purrr)
library(tibble)

calc_group <- function(idx){
  xmat <- MEs2[idx, , drop=FALSE]
  y <- md2$inflammation_score[idx]
  map_dfr(colnames(xmat), \(m){
    x <- xmat[, m]
    ok <- is.finite(x) & is.finite(y)
    n <- sum(ok)
    if(n < 3) return(tibble(module=m, n=n, cor=NA_real_, p=NA_real_))
    ct <- suppressWarnings(cor.test(x[ok], y[ok], method="spearman", exact=FALSE))
    tibble(module=m, n=n, cor=unname(ct$estimate), p=ct$p.value)
  })
}

res <- md2 %>%
  mutate(.cell = rownames(md2)) %>%
  group_by(T_subtype) %>%
  group_modify(\(df, key){
    idx <- match(df$.cell, rownames(MEs2))
    out <- calc_group(idx)
    out
  }) %>%
  ungroup() %>%
  group_by(T_subtype) %>%
  mutate(fdr = p.adjust(p, method="fdr")) %>%
  ungroup()


library(ggplot2)

res2 <- res %>%
  filter(module != "grey") %>%
  mutate(star = case_when(
    fdr < 0.001 ~ "***",
    fdr < 0.01 ~ "**",
    fdr < 0.05 ~ "*",
    TRUE ~ ""
  ))

p <- ggplot(res2, aes(x=T_subtype, y=module, fill=cor)) +
  geom_tile() +
  geom_text(aes(label=star), size=3) +
  theme_bw() +
  theme(axis.text.x = element_text(angle=45, hjust=1)) +
  labs(x="T_subtype", y="Module (hME)", fill="Spearman r",
       title="Module vs inflammation_score (only cells with score)")
setwd("data\\UC")  
p
ggsave("Module_inflammation_tile.pdf", p, width=10, height=8)
library(dplyr)

# 1) 模块基因表（包含每个基因的 module、kME 等）
mods <- GetModules(seurat_obj)
write.csv(mods, "treg_top_hubs.csv", row.names = FALSE)
# 看看模块名到底是 "treg5" 还是 "Treg5"
levels(mods$module)

# 2) 取 treg5 的基因，并按 kME 从大到小排
treg3_genes <- mods %>%
  filter(tolower(module) == "treg1") %>%   # 不区分大小写匹配
  arrange(desc(kME_Treg1)) %>%
  mutate(rank = row_number())

# 3) hub 基因：取前 N 个
n_hub <- 20
treg1_hubs <- treg3_genes %>% slice_head(n = n_hub)

treg1_hubs %>% select(gene_name, module, module, rank) %>% head(20)

# 导出
write.csv(treg5_hubs, "treg5_top_hubs.csv", row.names = FALSE)


####判断其和疾病的相关性####


seurat_obj1 <- ModuleTraitCorrelation(
  seurat_obj,
  traits = "disease",
  group.by = "T_subtype"
)
# 转 factor
seurat_obj$disease <- as.factor(seurat_obj$disease)

# 再跑
seurat_obj1 <- ModuleTraitCorrelation(
  seurat_obj,
  traits = "disease",
  group.by = "T_subtype"
)
MEs <- GetMEs(seurat_obj1, harmonized = TRUE)  # 建议用 harmonized
md  <- seurat_obj@meta.data

common <- intersect(rownames(md), rownames(MEs))
md2 <- md[common, , drop=FALSE]
MEs2 <- MEs[common, , drop=FALSE]

df <- cbind(md2, MEs2)
mod_cols <- colnames(MEs2)
head(mod_cols, 10)
library(ggplot2)

# 找到对应列名（按你实际列名调整）
treg1_col <- grep("Treg1", mod_cols, ignore.case = TRUE, value = TRUE)[1]

p_vln <- ggplot(df, aes(x = disease, y = .data[[treg1_col]])) +
  geom_violin(trim=FALSE) +
  geom_boxplot(width=0.15, outlier.size=0.2) +
  facet_wrap(~T_subtype, scales="free_y") +
  theme_bw() +
  theme(axis.text.x = element_text(angle=45, hjust=1)) +
  labs(x="Disease", y=treg1_col, title="Module activity by disease (facet by T_subtype)")

p_vln
ggsave("disease_Treg1_module_violin.pdf", p_vln, width=12, height=7)
library(dplyr)
library(tidyr)
library(purrr)
library(tibble)

# 确保是二分类（多分类可以先选两组，或做ANOVA/Kruskal）
levels(df$disease)

# 你想比较的两组（按你的实际标签改）
g1 <- levels(df$disease)[1]
g2 <- levels(df$disease)[2]

calc_effect <- function(d, m){
  x <- d[[m]]
  grp <- d$disease
  ok <- is.finite(x) & !is.na(grp)
  d <- d[ok, , drop=FALSE]
  if(length(unique(d$disease)) < 2) return(tibble(module=m, n=0, effect=NA_real_, p=NA_real_))
  
  # 用中位数差作为效应量（g2 - g1）
  med1 <- median(d[[m]][d$disease == g1])
  med2 <- median(d[[m]][d$disease == g2])
  eff  <- med2 - med1
  
  p <- suppressWarnings(wilcox.test(d[[m]] ~ d$disease)$p.value)
  tibble(module=m, n=nrow(d), effect=eff, p=p)
}

res <- df %>%
  group_by(T_subtype) %>%
  group_modify(~map_dfr(mod_cols, \(m) calc_effect(.x, m))) %>%
  ungroup() %>%
  group_by(T_subtype) %>%
  mutate(fdr = p.adjust(p, method="fdr"),
         star = case_when(fdr < 0.001 ~ "***",
                          fdr < 0.01 ~ "**",
                          fdr < 0.05 ~ "*",
                          TRUE ~ "")) %>%
  ungroup()

p_tile <- ggplot(res, aes(x=T_subtype, y=module, fill=effect)) +
  geom_tile() +
  geom_text(aes(label=star), size=3) +
  theme_bw() +
  theme(axis.text.x = element_text(angle=45, hjust=1)) +
  labs(x="T_subtype", y="Module (hME)", fill=paste0("Median(", g2, ") - Median(", g1, ")"),
       title="Disease effect on module activity (per T_subtype)")

p_tile
ggsave("disease_module_effect_tile.pdf", p_tile, width=12, height=0.25*length(mod_cols)+3)
library(WGCNA)

library(Seurat)
library(ggplot2)

# 1) 取 CD8 Trm
cd8_trm <- subset(seurat_obj, subset = T_subtype == "Trm_CD8")
cd8_trm_WITH <- subset(cd8_trm, subset = has_inflam_score == "yes")
# 2) 取 KLRG1 表达（默认用 RNA assay 的 data 槽：log-normalized）
DefaultAssay(cd8_trm_WITH) <- "RNA"
gene_sclect <- 'ITGA1'
df <- FetchData(
  cd8_trm_WITH,
  vars = c('ZEB2', "inflammation_score", "orig.ident", "gsm", "nCount_RNA", "percent.mt")
)

# 3) 看看炎症分数是否在细胞层面有变化（如果几乎全一样，就不要只看细胞层面）
summary(df$inflammation_score)
length(unique(df$inflammation_score))

# 4) 线性回归
fit <- lm(ZEB2 ~ inflammation_score, data = df)
summary(fit)

# 5) 相关（Pearson 线性相关；Spearman 单调相关更稳健）
cor.test(df$ZEB2, df$inflammation_score, method = "pearson")
cor.test(df$ZEB2, df$inflammation_score, method = "spearman")

# 6) 可视化：散点 + 线性拟合
ggplot(df, aes(x = inflammation_score, y = ZEB2)) +
  geom_point(size = 0.6, alpha = 0.25) +
  geom_smooth(method = "lm", se = TRUE) +
  theme_classic() +
  labs(x = "Inflammation score", y = "KLRG1 (log-normalized)")


library(WGCNA)
library(hdWGCNA)
library(tidyverse)
library(cowplot)
library(patchwork)
library(qs)
library(Seurat)
# 基础包
library(Seurat)
library(dplyr)
library(Matrix)

# hdWGCNA
library(hdWGCNA)
library(WGCNA)
####cd4组合####
# 建议：WGCNA 多线程
enableWGCNAThreads()
set.seed(1234)
IBD_T_CELL_remove_dcs2 <- readRDS("data/IBD_T_CELL_remove_dcs2.rds")
seurat_obj <- subset(IBD_T_CELL_remove_dcs2,disease=='UC')

# 确认meta列存在
seurat_obj <- SetupForWGCNA(
  seurat_obj,
  gene_select = "fraction", # the gene selection approach
  fraction = 0.05, # fraction of cells that a gene needs to be expressed in order to be included
  wgcna_name = "tutorial" # the name of the hdWGCNA experiment
)
# construct metacells  in each group
seurat_obj <- MetacellsByGroups(
  seurat_obj = seurat_obj,
  group.by = c("T_subtype", "orig.ident"), # specify the columns in seurat_obj@meta.data to group by
  reduction = 'harmony', # select the dimensionality reduction to perform KNN on
  k = 25, # nearest-neighbors parameter
  max_shared = 10, # maximum number of shared cells between two metacells
  ident.group = 'T_subtype' # set the Idents of the metacell seurat object
)

# normalize metacell expression matrix:
seurat_obj <- NormalizeMetacells(seurat_obj)
seurat_obj <- SetDatExpr(
  seurat_obj,
  group_name = c( 'TFH', 'TH17', 'Treg'), # the name of the group of interest in the group.by column
  group.by='T_subtype', # the metadata column containing the cell type info. This same column should have also been used in MetacellsByGroups
  assay = 'RNA', # using RNA assay
  layer = 'data' # using normalized data
)
library(future)
plan(sequential)
library(WGCNA)
disableWGCNAThreads()
# 或者
allowWGCNAThreads(nThreads = 2)

# Test different soft powers:
seurat_obj <- TestSoftPowers(
  seurat_obj,
  networkType = 'signed' # you can also use "unsigned" or "signed hybrid"
)

# plot the results:
plot_list <- PlotSoftPowers(seurat_obj)
getwd()
file.access(getwd(), 2)
# assemble with patchwork
wrap_plots(plot_list, ncol=2)
power_table <- GetPowerTable(seurat_obj)
head(power_table)
tom_dir <- file.path(tempdir(), "TOM")
dir.create(tom_dir, recursive = TRUE, showWarnings = FALSE)

seurat_obj <- ConstructNetwork(
  seurat_obj,
  soft_power = 8,
  tom_outdir = tom_dir,
  tom_name = "cd4_cells",
  overwrite_tom = TRUE
)
seurat_obj <- ConstructNetwork(
  seurat_obj,
  soft_power = 8, # 自定义了4，如果是自动选择的话可能会是3
  tom_outdir = "TOM",
  tom_name = 'cd4_cells', # name of the topoligical overlap matrix written to disk
  overwrite_tom = TRUE # 允许覆盖已存在的同名文件
)
PlotDendrogram(seurat_obj, main='cd4_cells hdWGCNA Dendrogram')
TOM <- GetTOM(seurat_obj)
seurat_obj <- ModuleEigengenes(
  seurat_obj,
  group.by.vars="orig.ident"
)
# harmonized module eigengenes:
hMEs <- GetMEs(seurat_obj)

# module eigengenes:
MEs <- GetMEs(seurat_obj, harmonized=FALSE)
# compute eigengene-based connectivity (kME):
seurat_obj <- ModuleConnectivity(
  seurat_obj,
  group.by = 'T_subtype', group_name = 'cd4_cells'
)
# rename the modules
seurat_obj <- ResetModuleNames(
  seurat_obj,
  new_name = "cd4_cells_"
)
# plot genes ranked by kME for each module
p <- PlotKMEs(seurat_obj, ncol=5)

p
# get the module assignment table:
modules <- GetModules(seurat_obj) %>% subset(module != 'grey')

# show the first 6 columns:
head(modules[,1:6])
# get hub genes
hub_df <- GetHubGenes(seurat_obj, n_hubs = 10)

head(hub_df)
# compute gene scoring for the top 25 hub genes by kME for each module
# with UCell method
library(UCell)
seurat_obj <- ModuleExprScore(
  seurat_obj,
  n_genes = 25,
  method='UCell'
)
saveRDS(seurat_obj, file='hdWGCNA_object_cd4_cells.rds')

#seurat_obj <- hdWGCNA_object_TREG

####判断其后炎症的相关性####
seurat_obj$inflammation_score <- seurat_obj$characteristics_ch1.14

# 检查是否有对不上的
sum(is.na(seurat_obj$inflammation_score))
# 如 NA 很多，说明两个对象的 cell barcode 不一致，需要先统一/映射

# 确保是数值型

## 1) 计算模块-表型相关（按细胞亚群分别算；你这里用 T_subtype，能单独得到 Treg 的相关矩阵）
seurat_obj_meta <- seurat_obj@meta.data%>%
  mutate(
    inflammation_score = as.character(inflammation_score),
    # 如果里面可能有 "score: 3.2" 这种，建议用 parse_number；纯数字用 as.numeric 也行
    inflammation_score = readr::parse_number(inflammation_score)
    # inflammation_score = as.numeric(inflammation_score)  # 纯数字字符串时也可
  ) 
seurat_obj@meta.data <- seurat_obj_meta
## 0) 保证 traits 是 numeric（NA 保留没关系）
x <- as.character(seurat_obj$inflammation_score)
x[x %in% c("NA","NaN","Inf","-Inf","", "NULL")] <- NA
x <- gsub(",", "", x)
seurat_obj$inflammation_score <- as.numeric(x)
seurat_obj$inflammation_score[!is.finite(seurat_obj$inflammation_score)] <- NA

## 1) group.by 必须是 factor（建议显式转一下）
seurat_obj$T_subtype <- as.factor(seurat_obj$T_subtype)

## 2) 做一个“是否有评分”的筛选列
seurat_obj$has_inflam_score <- factor(
  ifelse(is.na(seurat_obj$inflammation_score), "no", "yes"),
  levels = c("no", "yes")
)
# 确保是 numeric
seurat_obj$inflammation_score <- as.numeric(as.character(seurat_obj$inflammation_score))
seurat_obj$inflammation_score[!is.finite(seurat_obj$inflammation_score)] <- NA

# 用函数内部 subset（不要先 subset()）
seurat_obj <- ModuleTraitCorrelation(
  seurat_obj,
  traits = "inflammation_score",
  group.by = "T_subtype",
  subset_by = "has_inflam_score",
  subset_groups = "yes"
)










# 取 harmonized MEs（metacell层面通常用这个）
MEs <- GetMEs(seurat_obj, harmonized = TRUE)

# MEs 的行名应该是细胞/ metacell 的名字；把它拼到 meta.data
md <- seurat_obj@meta.data

# 对齐（取交集，避免不匹配）
common <- intersect(rownames(md), rownames(MEs))
md2 <- md[common, , drop=FALSE]
MEs2 <- MEs[common, , drop=FALSE]

# 只保留有炎症评分的
md2$inflammation_score <- as.numeric(as.character(md2$inflammation_score))
keep <- is.finite(md2$inflammation_score)
md2 <- md2[keep, , drop=FALSE]
MEs2 <- MEs2[keep, , drop=FALSE]


library(dplyr)
library(purrr)
library(tibble)

calc_group <- function(idx){
  xmat <- MEs2[idx, , drop=FALSE]
  y <- md2$inflammation_score[idx]
  map_dfr(colnames(xmat), \(m){
    x <- xmat[, m]
    ok <- is.finite(x) & is.finite(y)
    n <- sum(ok)
    if(n < 3) return(tibble(module=m, n=n, cor=NA_real_, p=NA_real_))
    ct <- suppressWarnings(cor.test(x[ok], y[ok], method="spearman", exact=FALSE))
    tibble(module=m, n=n, cor=unname(ct$estimate), p=ct$p.value)
  })
}

res <- md2 %>%
  mutate(.cell = rownames(md2)) %>%
  group_by(T_subtype) %>%
  group_modify(\(df, key){
    idx <- match(df$.cell, rownames(MEs2))
    out <- calc_group(idx)
    out
  }) %>%
  ungroup() %>%
  group_by(T_subtype) %>%
  mutate(fdr = p.adjust(p, method="fdr")) %>%
  ungroup()


library(ggplot2)

res2 <- res %>%
  filter(module != "grey") %>%
  mutate(star = case_when(
    fdr < 0.001 ~ "***",
    fdr < 0.01 ~ "**",
    fdr < 0.05 ~ "*",
    TRUE ~ ""
  ))

p <- ggplot(res2, aes(x=T_subtype, y=module, fill=cor)) +
  geom_tile() +
  geom_text(aes(label=star), size=3) +
  theme_bw() +
  theme(axis.text.x = element_text(angle=45, hjust=1)) +
  labs(x="T_subtype", y="Module (hME)", fill="Spearman r",
       title="Module vs inflammation_score (only cells with score)")
setwd("data\\UC")  
p
ggsave("Module_inflammation_tile.pdf", p, width=10, height=8)
library(dplyr)

# 1) 模块基因表（包含每个基因的 module、kME 等）
mods <- GetModules(seurat_obj)
write.csv(mods, "treg_top_hubs.csv", row.names = FALSE)
# 看看模块名到底是 "treg5" 还是 "Treg5"
levels(mods$module)

# 2) 取 treg5 的基因，并按 kME 从大到小排
treg3_genes <- mods %>%
  filter(tolower(module) == "treg1") %>%   # 不区分大小写匹配
  arrange(desc(kME_Treg1)) %>%
  mutate(rank = row_number())

# 3) hub 基因：取前 N 个
n_hub <- 20
treg1_hubs <- treg3_genes %>% slice_head(n = n_hub)

treg1_hubs %>% select(gene_name, module, module, rank) %>% head(20)

# 导出
write.csv(treg5_hubs, "treg5_top_hubs.csv", row.names = FALSE)


####判断其和疾病的相关性####


seurat_obj1 <- ModuleTraitCorrelation(
  seurat_obj,
  traits = "disease",
  group.by = "T_subtype"
)
# 转 factor
seurat_obj$disease <- as.factor(seurat_obj$disease)

# 再跑
seurat_obj1 <- ModuleTraitCorrelation(
  seurat_obj,
  traits = "disease",
  group.by = "T_subtype"
)
MEs <- GetMEs(seurat_obj1, harmonized = TRUE)  # 建议用 harmonized
md  <- seurat_obj@meta.data

common <- intersect(rownames(md), rownames(MEs))
md2 <- md[common, , drop=FALSE]
MEs2 <- MEs[common, , drop=FALSE]

df <- cbind(md2, MEs2)
mod_cols <- colnames(MEs2)
head(mod_cols, 10)
library(ggplot2)

# 找到对应列名（按你实际列名调整）
treg1_col <- grep("Treg1", mod_cols, ignore.case = TRUE, value = TRUE)[1]

p_vln <- ggplot(df, aes(x = disease, y = .data[[treg1_col]])) +
  geom_violin(trim=FALSE) +
  geom_boxplot(width=0.15, outlier.size=0.2) +
  facet_wrap(~T_subtype, scales="free_y") +
  theme_bw() +
  theme(axis.text.x = element_text(angle=45, hjust=1)) +
  labs(x="Disease", y=treg1_col, title="Module activity by disease (facet by T_subtype)")

p_vln
ggsave("disease_Treg1_module_violin.pdf", p_vln, width=12, height=7)
library(dplyr)
library(tidyr)
library(purrr)
library(tibble)

# 确保是二分类（多分类可以先选两组，或做ANOVA/Kruskal）
levels(df$disease)

# 你想比较的两组（按你的实际标签改）
g1 <- levels(df$disease)[1]
g2 <- levels(df$disease)[2]

calc_effect <- function(d, m){
  x <- d[[m]]
  grp <- d$disease
  ok <- is.finite(x) & !is.na(grp)
  d <- d[ok, , drop=FALSE]
  if(length(unique(d$disease)) < 2) return(tibble(module=m, n=0, effect=NA_real_, p=NA_real_))
  
  # 用中位数差作为效应量（g2 - g1）
  med1 <- median(d[[m]][d$disease == g1])
  med2 <- median(d[[m]][d$disease == g2])
  eff  <- med2 - med1
  
  p <- suppressWarnings(wilcox.test(d[[m]] ~ d$disease)$p.value)
  tibble(module=m, n=nrow(d), effect=eff, p=p)
}

res <- df %>%
  group_by(T_subtype) %>%
  group_modify(~map_dfr(mod_cols, \(m) calc_effect(.x, m))) %>%
  ungroup() %>%
  group_by(T_subtype) %>%
  mutate(fdr = p.adjust(p, method="fdr"),
         star = case_when(fdr < 0.001 ~ "***",
                          fdr < 0.01 ~ "**",
                          fdr < 0.05 ~ "*",
                          TRUE ~ "")) %>%
  ungroup()

p_tile <- ggplot(res, aes(x=T_subtype, y=module, fill=effect)) +
  geom_tile() +
  geom_text(aes(label=star), size=3) +
  theme_bw() +
  theme(axis.text.x = element_text(angle=45, hjust=1)) +
  labs(x="T_subtype", y="Module (hME)", fill=paste0("Median(", g2, ") - Median(", g1, ")"),
       title="Disease effect on module activity (per T_subtype)")

p_tile
ggsave("disease_module_effect_tile.pdf", p_tile, width=12, height=0.25*length(mod_cols)+3)
library(WGCNA)

library(Seurat)
library(ggplot2)

# 1) 取 CD8 Trm
cd8_trm <- subset(seurat_obj, subset = T_subtype == "Trm_CD8")
cd8_trm_WITH <- subset(cd8_trm, subset = has_inflam_score == "yes")
# 2) 取 KLRG1 表达（默认用 RNA assay 的 data 槽：log-normalized）
DefaultAssay(cd8_trm_WITH) <- "RNA"
gene_sclect <- 'ITGA1'
df <- FetchData(
  cd8_trm_WITH,
  vars = c('ZEB2', "inflammation_score", "orig.ident", "gsm", "nCount_RNA", "percent.mt")
)

# 3) 看看炎症分数是否在细胞层面有变化（如果几乎全一样，就不要只看细胞层面）
summary(df$inflammation_score)
length(unique(df$inflammation_score))

# 4) 线性回归
fit <- lm(ZEB2 ~ inflammation_score, data = df)
summary(fit)

# 5) 相关（Pearson 线性相关；Spearman 单调相关更稳健）
cor.test(df$ZEB2, df$inflammation_score, method = "pearson")
cor.test(df$ZEB2, df$inflammation_score, method = "spearman")

# 6) 可视化：散点 + 线性拟合
ggplot(df, aes(x = inflammation_score, y = ZEB2)) +
  geom_point(size = 0.6, alpha = 0.25) +
  geom_smooth(method = "lm", se = TRUE) +
  theme_classic() +
  labs(x = "Inflammation score", y = "KLRG1 (log-normalized)")

