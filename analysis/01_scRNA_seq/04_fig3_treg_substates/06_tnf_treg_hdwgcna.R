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
  group.by = c("treg_sub"), # specify the columns in seurat_obj@meta.data to group by
  reduction = 'harmony', # select the dimensionality reduction to perform KNN on
  k = 15, # nearest-neighbors parameter
  max_shared = 20, # maximum number of shared cells between two metacells
  ident.group = 'treg_sub' # set the Idents of the metacell seurat object
)
unique(seurat_obj$treg_sub)
# normalize metacell expression matrix:
seurat_obj <- NormalizeMetacells(seurat_obj)
seurat_obj <- SetDatExpr(
  seurat_obj,
  group_name = "TNF_TREG", # the name of the group of interest in the group.by column
  group.by='treg_sub', # the metadata column containing the cell type info. This same column should have also been used in MetacellsByGroups
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
  group.by = 'treg_sub', group_name = 'TNF_TREG'
)
# rename the modules
seurat_obj <- ResetModuleNames(
  seurat_obj,
  new_name = "TNF"
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
saveRDS(seurat_obj, file='hdWGCNA_object_TNF_TREG.rds')
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
seurat_obj$CELL_DISE <- paste0(seurat_obj$treg_sub,'_',seurat_obj$disease)
p <- DotPlot(seurat_obj, features=mods, group.by = 'CELL_DISE')
p
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
DotPlot(TREG_SUBSET,features = c('NFKB1','IFNG','IL10','TNFRSF1B','TNFRSF18','TGFB1'),group.by = 'treg_sub')
write.csv(modules,file = 'modules_TNF_TREG.CSV')
aaa <- FindMarkers(object = TREG_SUBSET,ident.1 = 'TNF_TREG')
aaa <- aaa %>% filter(p_val<0.05,abs(avg_log2FC)>1)

aaa <- FindAllMarkers(TREG_SUBSET,group.by = "treg_sub",only.pos = T,min.pct = 0.1,logfc.threshold = 0.5)
write.csv(aaa,file = 'TNF_TREG_marker.CSV')
