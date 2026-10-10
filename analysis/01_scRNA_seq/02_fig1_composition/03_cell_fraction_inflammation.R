library(dplyr)
library(ggplot2)
sce.all.filt <- IBD_T_TREG_SUB 
head(sce.all.filt@meta.data)
sce.all.filt$inflammation_score <- sce.all.filt$characteristics_ch1.14
analyze_subtype_vs_inflammation <- function(
    md,
    subtype,
    level_n = 4,
    method = "spearman"
) {
  
  df <- md %>%
    filter(!is.na(orig.ident), !is.na(T_subtype)) %>%
    mutate(
      inflammation_score = as.character(inflammation_score),
      # 如果里面可能有 "score: 3.2" 这种，建议用 parse_number；纯数字用 as.numeric 也行
      inflammation_score = readr::parse_number(inflammation_score)
      # inflammation_score = as.numeric(inflammation_score)  # 纯数字字符串时也可
    ) %>%
    group_by(orig.ident) %>%
    summarise(
      n_cells = n(),
      n_subtype = sum(T_subtype == subtype),
      prop_subtype = n_subtype / n_cells,
      inflammation_score = first(na.omit(inflammation_score)),
      .groups = "drop"
    ) %>%
    filter(!is.na(inflammation_score))
  
  if (nrow(df) < 5) {
    warning(paste(subtype, ": too few patients"))
    return(NULL)
  }
  
  # 防止所有病人的 inflammation_score 一样（cor.test 也会失败/没意义）
  if (dplyr::n_distinct(df$inflammation_score) < 2) {
    warning(paste(subtype, ": inflammation_score has <2 distinct values"))
    return(NULL)
  }
  if (dplyr::n_distinct(df$prop_subtype) < 2) {
    warning(paste(subtype, ": prop_subtype has <2 distinct values"))
    return(NULL)
  }
  
  cor_res <- cor.test(df$inflammation_score, df$prop_subtype, method = method)
  
  p_label <- ifelse(
    cor_res$p.value < 0.05,
    paste0("ρ = ", round(cor_res$estimate, 3), "\nP < 0.05", "\nn = ", nrow(df)),
    paste0("ρ = ", round(cor_res$estimate, 3), "\nP = ", signif(cor_res$p.value, 3), "\nn = ", nrow(df))
  )
  
  p <- ggplot(df, aes(x = inflammation_score, y = prop_subtype)) +
    geom_point(aes(size = n_cells), alpha = 0.7) +
    geom_smooth(method = "lm", se = TRUE) +
    annotate("text", x = Inf, y = Inf, label = p_label, hjust = 1.1, vjust = 1.1, size = 4) +
    labs(
      title = subtype,
      x = "Inflammation score",
      y = paste0("Proportion of ", subtype),
      size = "Cells per patient"
    ) +
    theme_bw()
  
  list(
    data = df,
    plot = p,
    stat = data.frame(
      subtype = subtype,
      n_patient = nrow(df),
      rho = as.numeric(cor_res$estimate),
      p_value = cor_res$p.value
    )
  )
}

md <- sce.all.filt@meta.data
subtypes <- sort(unique(md$T_subtype))
subtypes
res_list <- lapply(subtypes, function(st) {
  analyze_subtype_vs_inflammation(md, st)
})

# 去掉 NULL
res_list <- res_list[!sapply(res_list, is.null)]
stat_df <- bind_rows(lapply(res_list, `[[`, "stat"))

stat_df %>%
  arrange(p_value)
stat_df <- stat_df %>%
  mutate(FDR = p.adjust(p_value, method = "fdr"))

library(ggplot2)

ggplot(stat_df, aes(x = reorder(subtype, rho), y = rho, fill = FDR < 0.05)) +
  geom_col() +
  coord_flip() +
  labs(
    x = "Cell subtype",
    y = "Spearman ρ (vs inflammation)",
    fill = "FDR < 0.05"
  ) +
  theme_bw()
# 确保 age 是数值
sce.all.filt$age <- as.numeric(as.character(sce.all.filt$age))

# 分组：>35 为 old，否则 young
sce.all.filt$age_group <- ifelse(sce.all.filt$age >= 40, "old", "young")


# 看一下分布
table(sce.all.filt$age_group, useNA = "ifany")
sce.all.filt$age_group_sex <- paste0(sce.all.filt$age_group,'_',sce.all.filt$sex1)
sce.all.filt$age_group_sex_dis <- paste0(sce.all.filt$age_group_sex,'_',sce.all.filt$disease)
cellRatioPlot(object = DSS_MOUSE_T,
              sample.name = "condition",
              celltype.name = "T_subtype")
FeaturePlot(DSS_MOUSE_T,features = 'ITGAE')
DimPlot(DSS_MOUSE_T,group.by = 'T_subtype')


## 3) 画热图（可选）
p <- PlotModuleTraitCorrelation(
  seurat_obj,
  label = "fdr",
  label_symbol = "stars",
  text_size = 2,
  text_digits = 2,
  plot_max = 0.5,
  combine = TRUE
)  # 官方示例 :contentReference[oaicite:3]{index=3}
p
library(dplyr)
library(ggplot2)
library(readr)
library(scales)

## 如果 stat_df 还没有 FDR，先补上
if (!"FDR" %in% colnames(stat_df)) {
  stat_df <- stat_df %>%
    mutate(FDR = p.adjust(p_value, method = "fdr"))
}

## 1. 提取 Treg 每个病人的比例
df_treg <- md %>%
  filter(!is.na(orig.ident), !is.na(T_subtype)) %>%
  mutate(
    inflammation_score = readr::parse_number(as.character(inflammation_score))
  ) %>%
  group_by(orig.ident) %>%
  summarise(
    n_cells = n(),
    n_Treg = sum(T_subtype == "Treg"),
    prop_Treg = n_Treg / n_cells,
    inflammation_score = first(na.omit(inflammation_score)),
    .groups = "drop"
  ) %>%
  filter(!is.na(inflammation_score))

## 2. 提取 Treg 相关性统计结果
stat_treg <- stat_df %>%
  filter(subtype == "Treg")

label_treg <- paste0(
  "Spearman \u03c1 = ", round(stat_treg$rho, 3), "\n",
  "P = ", formatC(stat_treg$p_value, format = "e", digits = 2), "\n",
  "FDR = ", formatC(stat_treg$FDR, format = "e", digits = 2), "\n",
  "n = ", stat_treg$n_patient
)
p_treg <- ggplot(df_treg, aes(x = inflammation_score, y = prop_Treg)) +
  geom_point(
    color = "#0072B2",
    alpha = 0.75,
    size = 2.5
  ) +
  geom_smooth(
    method = "lm",
    se = TRUE,
    color = "#D55E00",
    fill = "#D55E00",
    alpha = 0.18,
    linewidth = 0.9
  ) +
  annotate(
    "text",
    x = Inf,
    y = Inf,
    label = label_treg,
    hjust = 1.05,
    vjust = 1.15,
    size = 4
  ) +
  scale_y_continuous(labels = percent_format(accuracy = 1)) +
  labs(
    title = "Treg",
    x = "Inflammation score",
    y = "Proportion of Treg"
  ) +
  theme_classic(base_size = 13) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5, size = 15),
    axis.title = element_text(face = "bold"),
    axis.text = element_text(color = "black")
  )

p_treg
library(dplyr)
library(ggplot2)
library(readr)
library(scales)
library(ggpubr)  # 用于正态性检验

## 如果 stat_df 还没有 FDR，先补上
if (!"FDR" %in% colnames(stat_df)) {
  stat_df <- stat_df %>%
    mutate(FDR = p.adjust(p_value, method = "fdr"))
}

## 1. 提取 Treg 每个病人的比例
df_treg <- md %>%
  filter(!is.na(orig.ident), !is.na(T_subtype)) %>%
  mutate(
    inflammation_score = readr::parse_number(as.character(inflammation_score))
  ) %>%
  group_by(orig.ident) %>%
  summarise(
    n_cells = n(),
    n_Treg = sum(T_subtype == "Treg"),
    prop_Treg = n_Treg / n_cells,
    inflammation_score = first(na.omit(inflammation_score)),
    .groups = "drop"
  ) %>%
  filter(!is.na(inflammation_score))

## 2. 正态性检验
# 检查 prop_Treg 和 inflammation_score 的正态性
shapiro_test_prop <- shapiro.test(df_treg$prop_Treg)
shapiro_test_inflam <- shapiro.test(df_treg$inflammation_score)

# 检查残差的正态性（在拟合线性模型后）
lm_fit <- lm(prop_Treg ~ inflammation_score, data = df_treg)
shapiro_test_residuals <- shapiro.test(residuals(lm_fit))

# 创建正态性检验结果标签
normality_label <- paste0(
  "Normality Tests:\n",
  "prop_Treg Shapiro-Wilk: p = ", formatC(shapiro_test_prop$p.value, format = "e", digits = 2), "\n",
  "inflammation_score Shapiro-Wilk: p = ", formatC(shapiro_test_inflam$p.value, format = "e", digits = 2), "\n",
  "Residuals Shapiro-Wilk: p = ", formatC(shapiro_test_residuals$p.value, format = "e", digits = 2), "\n",
  "n = ", nrow(df_treg)
)

## 3. 提取 Treg 相关性统计结果
stat_treg <- stat_df %>%
  filter(subtype == "Treg")

label_treg <- paste0(
  "Spearman \u03c1 = ", round(stat_treg$rho, 3), "\n",
  "P = ", formatC(stat_treg$p_value, format = "e", digits = 2), "\n",
  "FDR = ", formatC(stat_treg$FDR, format = "e", digits = 2), "\n",
  "n = ", stat_treg$n_patient
)

## 4. 绘制散点图和回归线
p_treg <- ggplot(df_treg, aes(x = inflammation_score, y = prop_Treg)) +
  geom_point(
    color = "#0072B2",
    alpha = 0.75,
    size = 2.5
  ) +
  geom_smooth(
    method = "lm",
    se = TRUE,
    color = "#D55E00",
    fill = "#D55E00",
    alpha = 0.18,
    linewidth = 0.9
  ) +
  annotate(
    "text",
    x = Inf,
    y = Inf,
    label = normality_label,
    hjust = 1.05,
    vjust = 0.85,
    size = 4,
    color = "darkred"
  ) +
  scale_y_continuous(labels = percent_format(accuracy = 1)) +
  labs(
    title = "Treg",
    x = "Inflammation score",
    y = "Proportion of Treg"
  ) +
  theme_classic(base_size = 13) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5, size = 15),
    axis.title = element_text(face = "bold"),
    axis.text = element_text(color = "black")
  )
p_treg 
## 5. 添加Q-Q图（可选）
# 如果您想单独查看Q-Q图
qq_plot <- ggqqplot(df_treg$prop_Treg, title = "Q-Q Plot of Treg Proportion") +
  theme_classic(base_size = 13)

# 组合主图和Q-Q图
combined_plot <- ggarrange(p_treg, qq_plot, ncol = 1, nrow = 2, heights = c(2, 1))

# 保存图像
ggsave("treg_analysis_with_normality_test.pdf", combined_plot, 
       width = 183, height = 183, units = "mm", dpi = 300)
library(dplyr)
library(ggplot2)
library(readr)
library(scales)
library(ggExtra)  # 添加ggExtra包

## 如果 stat_df 还没有 FDR，先补上
if (!"FDR" %in% colnames(stat_df)) {
  stat_df <- stat_df %>%
    mutate(FDR = p.adjust(p_value, method = "fdr"))
}

## 1. 提取 Treg 每个病人的比例
df_treg <- md %>%
  filter(!is.na(orig.ident), !is.na(T_subtype)) %>%
  mutate(
    inflammation_score = readr::parse_number(as.character(inflammation_score))
  ) %>%
  group_by(orig.ident) %>%
  summarise(
    n_cells = n(),
    n_Treg = sum(T_subtype == "Treg"),
    prop_Treg = n_Treg / n_cells,
    inflammation_score = first(na.omit(inflammation_score)),
    .groups = "drop"
  ) %>%
  filter(!is.na(inflammation_score))

## 2. 提取 Treg 相关性统计结果
stat_treg <- stat_df %>%
  filter(subtype == "Treg")

label_treg <- paste0(
  "Spearman \u03c1 = ", round(stat_treg$rho, 3), "\n",
  "P = ", formatC(stat_treg$p_value, format = "e", digits = 2), "\n",
  "FDR = ", formatC(stat_treg$FDR, format = "e", digits = 2), "\n",
  "n = ", stat_treg$n_patient
)

## 3. 创建基础散点图
p_treg_base <- ggplot(df_treg, aes(x = inflammation_score, y = prop_Treg)) +
  geom_point(
    color = "#0072B2",
    alpha = 0.75,
    size = 2.5
  ) +
  geom_smooth(
    method = "lm",
    se = TRUE,
    color = "#D55E00",
    fill = "#D55E00",
    alpha = 0.18,
    linewidth = 0.9
  ) +
  scale_y_continuous(labels = percent_format(accuracy = 1)) +
  labs(
    title = "Treg",
    x = "Inflammation score",
    y = "Proportion of Treg"
  ) +
  theme_classic(base_size = 13) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5, size = 15),
    axis.title = element_text(face = "bold"),
    axis.text = element_text(color = "black")
  )

## 4. 添加边际直方图
p_treg <- ggMarginal(
  p_treg_base,
  type = "histogram",  # 可以选择 "density", "boxplot", "violin" 等
  margins = "both",    # "x", "y", 或 "both"
  size = 5,           # 直方图大小
  fill = "gray80",     # 直方图填充颜色
  color = "black",    # 直方图边框颜色
  xparams = list(bins = 20),  # x轴直方图参数
  yparams = list(bins = 15)   # y轴直方图参数
) +
  # 添加统计标签
  annotate(
    "text",
    x = Inf,
    y = Inf,
    label = label_treg,
    hjust = 1.05,
    vjust = 1.15,
    size = 4
  )

## 5. 显示图形
p_treg

## 6.
p_treg_custom <- ggMarginal(
  p_treg_base,
  type = "histogram",
  margins = "both",
  size = 5,
  fill = "gray80",
  color = "black",
  xparams = list(
    bins = 25,
    position = "identity",
    alpha = 0.7
  ),
  yparams = list(
    bins = 25,
    position = "identity",
    alpha = 0.7
  ),
  ggproto = list(
    # 自定义直方图外观
    geom_bar = function(data, params) {
      ggplot2::geom_bar(
        data = data,
        stat = "bin",
        position = params$position,
        alpha = params$alpha,
        fill = params$fill,
        color = params$color,
        width = 0.9
      )
    }
  )
)
p_treg_custom
