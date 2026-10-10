paste0(IBD_T_CELL_remove_dcs2$CELL_DISE)
IBD_T_CELL_remove_dcs2$CELL_DISE <- paste0(IBD_T_CELL_remove_dcs2$T_subtype,'_',IBD_T_CELL_remove_dcs2$disease)
DotPlot(IBD_T_CELL_remove_dcs2,features = c('BATF','FOS','IRF4'),group.by = 'CELL_DISE')+RotatedAxis()


library(Seurat)
library(dplyr)
library(tidyr)
library(purrr)
####病人层面####
# 1. 只取 Treg
treg_seurat <- subset(IBD_T_CELL_remove_dcs2, T_subtype == "TFH")

# 2. 取要分析的基因
genes_to_test <- 'BATF'

# 保险起见，只保留 Seurat 对象中真的存在的基因
genes_to_test <- intersect(genes_to_test, rownames(treg_seurat))


# 3. 提取表达矩阵（基因 x 细胞）
expr_mat <- GetAssayData(treg_seurat, slot = "data")[genes_to_test, , drop = FALSE]

# 转成 data.frame：行为细胞，列为基因
expr_df <- as.data.frame(t(as.matrix(expr_mat)))
expr_df$cell <- rownames(expr_df)

# 4. 提取 metadata
meta_df <- treg_seurat@meta.data %>%
  mutate(
    cell = rownames(.),
    inflammation_score_num = ifelse(
      disease == "HC",
      0,
      as.numeric(gsub("inflammation score: ", "", as.character(inflammation_score)))
    )
  ) %>%
  select(cell, orig.ident, disease, inflammation_score_num)

# 5. 合并表达和 metadata
cell_data <- left_join(meta_df, expr_df, by = "cell") %>%
  filter(!is.na(inflammation_score_num))

# 6. 按病人汇总：每个病人的每个基因平均表达
patient_expr <- cell_data %>%
  group_by(orig.ident) %>%
  summarise(
    inflammation_score = unique(inflammation_score_num)[1],
    disease = unique(disease)[1],
    across(all_of(genes_to_test), ~ mean(.x, na.rm = TRUE)),
    n_cells = n(),
    .groups = "drop"
  )

# 看一下病人层面数据
head(patient_expr)
# 7. 逐个基因做病人层面相关分析
cor_results <- map_dfr(genes_to_test, function(gene) {
  x <- patient_expr[[gene]]
  y <- patient_expr$inflammation_score
  
  # 去掉 NA
  keep <- !is.na(x) & !is.na(y)
  x <- x[keep]
  y <- y[keep]
  
  # 样本太少时跳过
  if (length(x) < 3 || length(unique(x)) < 2 || length(unique(y)) < 2) {
    return(data.frame(
      gene = gene,
      n_patients = length(x),
      spearman_rho = NA,
      spearman_p = NA,
      pearson_r = NA,
      pearson_p = NA
    ))
  }
  
  spearman_res <- suppressWarnings(cor.test(x, y, method = "spearman"))
  pearson_res  <- suppressWarnings(cor.test(x, y, method = "pearson"))
  
  data.frame(
    gene = gene,
    n_patients = length(x),
    spearman_rho = as.numeric(spearman_res$estimate),
    spearman_p = spearman_res$p.value,
    pearson_r = as.numeric(pearson_res$estimate),
    pearson_p = pearson_res$p.value
  )
})

# 8. 多重检验校正
cor_results <- cor_results %>%
  mutate(
    spearman_padj = p.adjust(spearman_p, method = "BH"),
    pearson_padj  = p.adjust(pearson_p, method = "BH")
  ) %>%
  arrange(desc(abs(spearman_rho)))

cor_results <- cor_results %>% filter(spearman_padj<0.05&pearson_padj<0.05)

library(ggplot2)

top_genes <- cor_results %>%
  arrange(desc(abs(spearman_rho))) %>%
  slice(1:6) %>%
  pull(gene)

plot_df <- patient_expr %>%
  select(orig.ident, inflammation_score, all_of('BATF')) %>%
  pivot_longer(cols = all_of('BATF'), names_to = "gene", values_to = "expr")

cor_annot <- cor_results %>%
  mutate(label = paste0("rho = ", round(spearman_rho, 2),
                        "\n p = ", signif(spearman_p, 3)))

# 2. 绘制散点 + 回归线 + facet_wrap
ggplot(plot_df, aes(x = inflammation_score, y = expr)) +
  geom_point(size = 2.5, color = "#4C8BF5") +
  geom_smooth(method = "lm", se = TRUE, color = "red", fill = "pink") +
  facet_wrap(~ gene, scales = "free_y") +
  theme_classic() +
  labs(
    x = "Inflammation Score",
    y = "Mean Expression per Patient",
    title = "Patient-level correlation between BATF of TFH and inflammation score"
  ) +
  # 3. 添加 Spearman rho 和 p 值注释
  geom_text(
    data = cor_annot,
    aes(
      x = max(plot_df$inflammation_score, na.rm = TRUE) *0.9,  # 调整位置
      y = max(plot_df$expr, na.rm = TRUE) *1,
      label = label
    ),
    inherit.aes = FALSE,
    hjust = 0,
    size = 4,
    fontface = "bold"
  )


####细胞层面####
library(Seurat)
library(dplyr)
library(purrr)
library(tidyr)
library(ggplot2)

# 1. 只取 UC 的 Treg
treg_seurat <- subset(IBD_T_CELL_remove_dcs2, T_subtype == "Treg" & disease == "UC")

# 2. 取要分析的基因
genes_to_test <- treg_wgcna$gene_name
genes_to_test <- intersect(genes_to_test, rownames(treg_seurat))

# 3. 提取表达矩阵
expr_mat <- GetAssayData(treg_seurat, slot = "data")[genes_to_test, , drop = FALSE]

# 转成 data.frame：行为细胞，列为基因
expr_df <- as.data.frame(t(as.matrix(expr_mat)))
expr_df$cell <- rownames(expr_df)

# 4. 提取 metadata
meta_df <- treg_seurat@meta.data %>%
  mutate(
    cell = rownames(.),
    inflammation_score_num = as.numeric(gsub("inflammation score: ", "", as.character(inflammation_score)))
  ) %>%
  select(cell, orig.ident, disease, inflammation_score_num)

# 5. 合并表达和 metadata
cell_data <- left_join(meta_df, expr_df, by = "cell") %>%
  filter(!is.na(inflammation_score_num))

head(cell_data)
cor_results <- map_dfr(genes_to_test, function(gene) {
  x <- cell_data[[gene]]
  y <- cell_data$inflammation_score_num
  
  # 去掉 NA，并去掉该基因表达为 0 的细胞
  keep <- !is.na(x) & !is.na(y) & x > 0
  x <- x[keep]
  y <- y[keep]
  
  if (length(x) < 3 || length(unique(x)) < 2 || length(unique(y)) < 2) {
    return(data.frame(
      gene = gene,
      n_cells = length(x),
      spearman_rho = NA,
      spearman_p = NA,
      pearson_r = NA,
      pearson_p = NA
    ))
  }
  
  spearman_res <- suppressWarnings(cor.test(x, y, method = "spearman"))
  pearson_res  <- suppressWarnings(cor.test(x, y, method = "pearson"))
  
  data.frame(
    gene = gene,
    n_cells = length(x),
    spearman_rho = as.numeric(spearman_res$estimate),
    spearman_p = spearman_res$p.value,
    pearson_r = as.numeric(pearson_res$estimate),
    pearson_p = pearson_res$p.value
  )
})

cor_results <- cor_results %>%
  mutate(
    spearman_padj = p.adjust(spearman_p, method = "BH"),
    pearson_padj  = p.adjust(pearson_p, method = "BH")
  ) %>%
  arrange(desc(abs(spearman_rho)))

cor_results
top_genes <- cor_results %>%
  arrange(desc(abs(spearman_rho))) %>%
  slice(1:6) %>%
  pull(gene)

plot_df <- cell_data %>%
  select(inflammation_score_num, all_of(top_genes)) %>%
  pivot_longer(cols = all_of(top_genes), names_to = "gene", values_to = "expr") %>%
  filter(expr > 0)
ggplot(plot_df, aes(x = inflammation_score_num, y = expr)) +
  geom_jitter(width = 0.1, height = 0, alpha = 0.25, size = 0.8, color = "#4C8BF5") +
  geom_smooth(method = "lm", se = TRUE, color = "red", fill = "pink") +
  facet_wrap(~ gene, scales = "free_y") +
  theme_classic() +
  labs(
    x = "Inflammation Score",
    y = "Expression per Cell (>0 only)",
    title = "Cell-level correlation between Treg genes and inflammation score"
  )

####细胞比例层面####
library(dplyr)filter(!is.na(inflammation_score_num))
library(ggplot2)

SUBSET <- subset(IBD_T_CELL_remove_dcs2, subset = !is.na(inflammation_score))
meta <- IBD_T_CELL_remove_dcs2@meta.data %>%
  mutate(
    cell = rownames(.),
    inflammation_score_num = ifelse(
      disease == "HC",
      0,
      as.numeric(gsub("inflammation score: ", "", as.character(inflammation_score)))
    )
  )

# 每个病人的总T细胞数、Treg数、Th17数
patient_prop <- meta %>%
  group_by(orig.ident) %>%
  summarise(
    disease = unique(disease)[1],
    inflammation_score = unique(inflammation_score_num)[1],
    total_T = n(),
    n_Treg = sum(T_subtype == "Treg", na.rm = TRUE),
    n_Th17 = sum(T_subtype == "TH17", na.rm = TRUE),
    Treg_prop = n_Treg / total_T,
    Th17_prop = n_Th17 / total_T,
    Treg_Th17_ratio = n_Treg / ifelse(n_Th17 == 0, NA, n_Th17),
    Treg_fraction_in_TregTh17 = n_Treg / ifelse((n_Treg + n_Th17) == 0, NA, (n_Treg + n_Th17)),
    .groups = "drop"
  )

patient_prop
cor.test(patient_prop$inflammation_score, patient_prop$Treg_prop, method = "spearman")
cor.test(patient_prop$inflammation_score, patient_prop$Treg_prop, method = "pearson")
cor.test(patient_prop$inflammation_score, patient_prop$Th17_prop, method = "spearman")
cor.test(patient_prop$inflammation_score, patient_prop$Th17_prop, method = "pearson")
ratio_df <- patient_prop %>% filter(!is.na(Treg_Th17_ratio))

cor.test(ratio_df$inflammation_score, ratio_df$Treg_Th17_ratio, method = "spearman")
cor.test(ratio_df$inflammation_score, ratio_df$Treg_Th17_ratio, method = "pearson")



cor_treg  <- cor.test(patient_prop$inflammation_score, patient_prop$Treg_prop, method = "pearson")
cor_th17  <- cor.test(patient_prop$inflammation_score, patient_prop$Th17_prop, method = "pearson")
cor_ratio <- cor.test(ratio_df$inflammation_score, ratio_df$Treg_Th17_ratio, method = "pearson")
rho_treg <- round(as.numeric(cor_treg$estimate), 3)
p_treg   <- signif(cor_treg$p.value, 3)

rho_th17 <- round(as.numeric(cor_th17$estimate), 3)
p_th17   <- signif(cor_th17$p.value, 3)

rho_ratio <- round(as.numeric(cor_ratio$estimate), 3)
p_ratio   <- signif(cor_ratio$p.value, 3)


ggplot(patient_prop, aes(x = inflammation_score, y = Treg_prop)) +
  geom_point(size = 3, color = "#4C8BF5") +
  geom_smooth(method = "lm", se = TRUE, color = "red", fill = "pink") +
  annotate(
    "text",
    x = Inf, y = Inf,
    label = paste0("Spearman rho = ", rho_treg, "\nP = ", p_treg),
    hjust = 1.1, vjust = 1.5, size = 5
  ) +
  theme_classic() +
  labs(
    title = "Treg proportion vs inflammation score",
    x = "Inflammation score",
    y = "Treg proportion"
  )
ggplot(patient_prop, aes(x = inflammation_score, y = Th17_prop)) +
  geom_point(size = 3, color = "#4C8BF5") +
  geom_smooth(method = "lm", se = TRUE, color = "red", fill = "pink") +
  annotate(
    "text",
    x = Inf, y = Inf,
    label = paste0("Spearman rho = ", rho_th17, "\nP = ", p_th17),
    hjust = 1.1, vjust = 1.5, size = 5
  ) +
  theme_classic() +
  labs(
    title = "Th17 proportion vs inflammation score",
    x = "Inflammation score",
    y = "Th17 proportion"
  )
ggplot(ratio_df, aes(x = inflammation_score, y = Treg_Th17_ratio)) +
  geom_point(size = 3, color = "#4C8BF5") +
  geom_smooth(method = "lm", se = TRUE, color = "red", fill = "pink") +
  annotate(
    "text",
    x = Inf, y = Inf,
    label = paste0("Spearman rho = ", rho_ratio, "\nP = ", p_ratio),
    hjust = 1.1, vjust = 1.5, size = 5
  ) +
  theme_classic() +
  labs(
    title = "Treg/Th17 ratio vs inflammation score",
    x = "Inflammation score",
    y = "Treg / Th17 ratio"
  )

