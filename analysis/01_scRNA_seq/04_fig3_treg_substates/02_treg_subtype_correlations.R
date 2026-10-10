library(Seurat)
library(dplyr)
library(ggplot2)
library(ggrepel)
table(IBD_T_TREG_SUB$T_subtype)
# 1. 提取 Treg
treg_cells <- subset(
  IBD_T_TREG_SUB,
  subset = T_subtype == "Treg"
)

# 检查 HC / UC 中 Treg 数量
table(treg_cells$disease)

# 2. 设置分组
Idents(treg_cells) <- "disease"

# 3. 差异分析：UC vs HC
deg_treg <- FindMarkers(
  object = treg_cells,
  ident.1 = "UC",
  ident.2 = "HC",
  assay = "RNA",
  slot = "data",
  
  logfc.threshold = 0,   # 不在分析前按FC筛基因，完整保留
  min.pct = 0.1,         # 任一组至少10%细胞表达
  min.diff.pct = -Inf,   # 不限制两组表达比例差异
  
  test.use = "wilcox",
  only.pos = FALSE,
  max.cells.per.ident = Inf,
  verbose = TRUE
)

deg_treg$gene <- rownames(deg_treg)

# 兼容不同 Seurat 版本
fc_col <- if ("avg_log2FC" %in% colnames(deg_treg)) {
  "avg_log2FC"
} else {
  "avg_logFC"
}

# 查看 BATF 结果
deg_treg[deg_treg$gene == "BATF", ]

# 4. 阈值：放宽至 FDR < 0.1
logfc_cutoff <- 0.5
padj_cutoff <- 0.1

deg_treg$p_val_adj[deg_treg$p_val_adj == 0] <- 1e-300

deg_treg$group <- "Not significant"

deg_treg$group[
  deg_treg[[fc_col]] > logfc_cutoff &
    deg_treg$p_val_adj < padj_cutoff
] <- "Up in UC"

deg_treg$group[
  deg_treg[[fc_col]] < -logfc_cutoff &
    deg_treg$p_val_adj < padj_cutoff
] <- "Up in HC"

# 取两侧显著基因用于标注
top_up <- deg_treg %>%
  filter(group == "Up in UC") %>%
  arrange(desc(.data[[fc_col]])) %>%
  slice_head(n = 8)

top_down <- deg_treg %>%
  filter(group == "Up in HC") %>%
  arrange(.data[[fc_col]]) %>%
  slice_head(n = 8)

# BATF 强制加入标注
label_genes <- deg_treg %>%
  filter(gene %in% c(top_up$gene, top_down$gene, "BATF"))

# 5. 火山图
p_volcano_treg <- ggplot(
  deg_treg,
  aes(
    x = .data[[fc_col]],
    y = -log10(p_val_adj),
    color = group
  )
) +
  geom_point(size = 1.25, alpha = 0.7) +
  scale_color_manual(
    values = c(
      "Up in UC" = "#D73027",
      "Up in HC" = "#4575B4",
      "Not significant" = "grey75"
    )
  ) +
  geom_vline(
    xintercept = c(-logfc_cutoff, logfc_cutoff),
    linetype = "dashed",
    linewidth = 0.4
  ) +
  geom_hline(
    yintercept = -log10(padj_cutoff),
    linetype = "dashed",
    linewidth = 0.4
  ) +
  geom_point(
    data = subset(deg_treg, gene == "BATF"),
    shape = 21,
    size = 4,
    fill = "yellow",
    color = "black",
    stroke = 1.1
  ) +
  geom_text_repel(
    data = label_genes,
    aes(label = gene),
    size = 4,
    fontface = "bold",
    max.overlaps = Inf
  ) +
  labs(
    title = "Treg differential expression: UC vs HC",
    subtitle = "T_subtype = Treg；阈值：|log2FC| > 0.5，FDR < 0.1",
    x = "log2 Fold Change (UC vs HC)",
    y = "-log10 adjusted P value",
    color = NULL
  ) +
  theme_classic(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5),
    legend.position = "top"
  )

p_volcano_treg



library(Seurat)
library(dplyr)
library(ggplot2)

# 1. 将字符串炎症评分转换为数值
treg_cells$inflammation_score_num <- as.numeric(
  sub(
    "^inflammation score:\\s*",
    "",
    as.character(treg_cells$inflammation_score)
  )
)

summary(treg_cells$inflammation_score_num)
table(is.na(treg_cells$inflammation_score_num))
sample_col <- "orig.ident"  # 改成你自己的样本编号列名

stopifnot(sample_col %in% colnames(treg_cells@meta.data))
# 2. 提取BATF表达、炎症评分、样本信息和疾病信息
dat <- FetchData(
  treg_cells,
  vars = c("BATF", "inflammation_score_num", sample_col, "disease")
)

colnames(dat)[colnames(dat) == sample_col] <- "sample_id"

# 3. 每个样本汇总为一个独立观测值
sample_df <- dat %>%
  group_by(sample_id, disease) %>%
  summarise(
    BATF_mean = mean(BATF, na.rm = TRUE),
    inflammation_score = median(inflammation_score_num, na.rm = TRUE),
    n_cells = n(),
    
    # 每个样本中有炎症评分的细胞数
    n_score_cells = sum(!is.na(inflammation_score_num)),
    
    # 每个样本中炎症评分的不同取值数
    n_score = n_distinct(inflammation_score_num[!is.na(inflammation_score_num)]),
    
    .groups = "drop"
  ) %>%
  filter(!is.na(inflammation_score))

# 检查：若炎症评分是样本级变量，n_score理论上应主要为1
sample_df %>% select(sample_id, disease, inflammation_score, n_score, n_cells)

# 4. 仅限UC样本，避免HC/UC分组本身造成伪相关
sample_uc <- sample_df %>%
  filter(disease == "UC")

# 推荐Spearman相关：更适合非正态、非线性单细胞表达数据
cor_uc <- cor.test(
  sample_uc$BATF_mean,
  sample_uc$inflammation_score,
  method = "spearman",
  exact = FALSE,
  alternative = "two.sided"
)

cor_uc
sample_uc <- sample_df %>%
  filter(disease == "UC"&inflammation_score>2)

cor_uc <- cor.test(
  sample_uc$BATF_mean,
  sample_uc$inflammation_score,
  method = "spearman",
  exact = FALSE
)

cor_uc
rho <- round(unname(cor_uc$estimate), 3)
p_value <- signif(cor_uc$p.value, 3)

ggplot(
  sample_uc,
  aes(x = inflammation_score, y = BATF_mean)
) +
  geom_point(size = 3) +
  geom_smooth(method = "lm", se = TRUE) +
  annotate(
    "text",
    x = Inf,
    y = Inf,
    hjust = 1.1,
    vjust = 1.5,
    label = paste0("Spearman rho = ", rho, "\nP = ", p_value),
    size = 5
  ) +
  theme_classic(base_size = 14) +
  labs(
    x = "Inflammation score",
    y = "Mean normalized BATF expression",
    title = "Correlation between BATF and inflammation score in UC-Treg"
  )










library(dplyr)
library(ggplot2)

# 提取每个Treg细胞的BATF表达与炎症评分
cell_df <- FetchData(
  treg_cells,
  vars = c("BATF", "inflammation_score_num")
) %>%
  filter(
    !is.na(inflammation_score_num)
  )

# 如果你还想只保留 BATF > 0 的细胞，取消下一行注释
cell_df <- cell_df %>% filter(BATF > 0&inflammation_score_num>0)

# 查看保留下来的细胞数
nrow(cell_df)

# 细胞水平 Spearman 相关性
cor_cell <- cor.test(
  cell_df$BATF,
  cell_df$inflammation_score_num,
  method = "spearman",
  exact = FALSE
)

cor_cell
rho <- round(unname(cor_cell$estimate), 3)
p_value <- signif(cor_cell$p.value, 3)

p_batf_cor <- ggplot(
  cell_df,
  aes(x = inflammation_score_num, y = BATF)
) +
  geom_point(
    size = 0.5,
    alpha = 0.25
  ) +
  geom_smooth(
    method = "lm",
    se = TRUE,
    linewidth = 1
  ) +
  annotate(
    "text",
    x = Inf,
    y = Inf,
    hjust = 1.1,
    vjust = 1.5,
    label = paste0(
      "Spearman rho = ", rho,
      "\nP = ", p_value,
      "\nCells = ", nrow(cell_df)
    ),
    size = 5
  ) +
  theme_classic(base_size = 14) +
  labs(
    title = "Cell-level correlation between BATF and inflammation score",
    x = "Inflammation score (> 0)",
    y = "BATF normalized expression"
  )

p_batf_cor

cluster_colors <- c(
  B= "#e9262a",
  'other' = "#BC5428FF",
  Epithelial = "#F876C2FF",
  glia = "#E54675FF",
  Mast = "#20A369FF",
  MYLOID= "#F4F071FF",
  neuron= "#10A4F3FF",
  endothelial_cells = "#9E180AFF",
  pericyte= "#612877",
  plasma = "#2278b4",
  Stromal = "#a6cde1",
  'T/NK' = "#fdce89"
)
cellRatioPlot(object = major_annotation_ibd_subtype,
                             sample.name = "dataset",
                             celltype.name = "cell_subtype",fill.col = cluster_colors)
