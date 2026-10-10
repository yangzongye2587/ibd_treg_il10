library(Seurat)
library(dplyr)
library(tidyr)
library(ggplot2)

obj <- major_annotation_ibd_subtype_tregback

celltype_col <- "cell_subtype"
group_col    <- "disease"
sample_col   <- "orig.ident"

meta_raw <- obj@meta.data

meta <- data.frame(
  sample   = as.character(meta_raw[[sample_col]]),
  disease  = as.character(meta_raw[[group_col]]),
  cell_subtype = as.character(meta_raw[[celltype_col]]),
  stringsAsFactors = FALSE
)

meta <- meta %>%
  dplyr::filter(!is.na(sample), !is.na(disease), !is.na(cell_subtype))

meta$disease <- factor(meta$disease, levels = c("HC", "UC"))

# 检查一个 sample 是否只对应一个 disease
sample_check <- meta %>%
  dplyr::distinct(sample, disease) %>%
  dplyr::count(sample, name = "n_group")

if (any(sample_check$n_group > 1)) {
  stop("同一个 sample 对应了多个 disease。请确认 orig.ident 是否是真正的样本/病人 ID。")
}

# 样本分组信息
sample_info <- meta %>%
  dplyr::distinct(sample, disease)

# 每个 sample 中每种 cell_subtype 的细胞数
count_df <- meta %>%
  dplyr::count(sample, cell_subtype, name = "n")

# 补齐缺失组合：某个 sample 没有某类细胞时，记为 0
sample_prop_df <- tidyr::expand_grid(
  sample = unique(meta$sample),
  cell_subtype = unique(meta$cell_subtype)
) %>%
  dplyr::left_join(count_df, by = c("sample", "cell_subtype")) %>%
  dplyr::mutate(n = tidyr::replace_na(n, 0)) %>%
  dplyr::left_join(sample_info, by = "sample") %>%
  dplyr::group_by(sample) %>%
  dplyr::mutate(
    total_cells = sum(n),
    sample_prop = n / total_cells
  ) %>%
  dplyr::ungroup()

# 先按 disease 对 sample-level 比例求平均
group_mean_df <- sample_prop_df %>%
  dplyr::group_by(cell_subtype, disease) %>%
  dplyr::summarise(
    mean_sample_prop = mean(sample_prop, na.rm = TRUE),
    .groups = "drop"
  )

# 在每个 cell_subtype 内，对 HC 和 UC 进行归一化
# 这样每个 cell_subtype 的柱子总高度都是 1
plot_df <- group_mean_df %>%
  dplyr::group_by(cell_subtype) %>%
  dplyr::mutate(
    rel_num = mean_sample_prop / sum(mean_sample_prop, na.rm = TRUE)
  ) %>%
  dplyr::ungroup()

# 按 UC 相对比例排序，可改成 HC 或字母顺序
celltype_order <- plot_df %>%
  dplyr::filter(disease == "UC") %>%
  dplyr::arrange(dplyr::desc(rel_num)) %>%
  dplyr::pull(cell_subtype)

plot_df$cell_subtype <- factor(plot_df$cell_subtype, levels = celltype_order)

p <- ggplot(plot_df, aes(x = cell_subtype, y = rel_num, fill = disease)) +
  geom_col(
    position = "stack",
    color = "black",
    linewidth = 0.35,
    width = 0.7
  ) +
  geom_text(
    aes(label = ifelse(rel_num > 0.05, paste0(round(rel_num * 100, 1), "%"), "")),
    position = position_stack(vjust = 0.5),
    size = 3,
    color = "black"
  ) +
  scale_fill_manual(values = c("HC" = "#96cada", "UC" = "#fdb26f")) +
  scale_y_continuous(
    limits = c(0, 1),
    breaks = seq(0, 1, 0.25),
    labels = function(x) paste0(x * 100, "%"),
    expand = c(0, 0)
  ) +
  labs(
    x = NULL,
    y = "Normalized cell ratio",
    fill = NULL,
    title = "Normalized disease composition within each cell subtype",
    subtitle = "Proportions were first calculated within each sample, averaged by disease group, then normalized within each cell subtype"
  ) +
  theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 15, hjust = 0.5),
    plot.subtitle = element_text(size = 10, hjust = 0.5, color = "grey35"),
    axis.text.x = element_text(angle = 45, hjust = 1, color = "black"),
    axis.text.y = element_text(color = "black"),
    axis.title.y = element_text(face = "bold"),
    legend.position = "right",
    legend.text = element_text(size = 10),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
    plot.margin = margin(10, 12, 10, 10)
  )

p
unique(major_annotation_ibd_subtype_tregback$sex1)
unique(major_annotation_ibd_subtype_tregback$age)
unique(major_annotation_ibd_subtype_tregback$tissue)

library(Seurat)
library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)

obj <- major_annotation_ibd_subtype_tregback

sample_col <- "orig.ident"
group_col  <- "disease"
sex_col    <- "sex1"
age_col    <- "age"
tissue_col <- "tissue"

meta_raw <- obj@meta.data

meta_cov <- data.frame(
  sample  = as.character(meta_raw[[sample_col]]),
  disease = as.character(meta_raw[[group_col]]),
  sex     = as.character(meta_raw[[sex_col]]),
  age     = as.character(meta_raw[[age_col]]),
  tissue  = as.character(meta_raw[[tissue_col]]),
  stringsAsFactors = FALSE
)

# =========================
# 1. 清洗 metadata
# =========================

meta_cov <- meta_cov %>%
  dplyr::mutate(
    disease = factor(disease, levels = c("HC", "UC")),
    sex = factor(sex, levels = c("Male", "Female")),
    age_num = suppressWarnings(as.numeric(age)),
    tissue_raw = trimws(tissue),
    tissue_lower = tolower(tissue_raw),
    tissue_clean = dplyr::case_when(
      is.na(tissue_raw) | tissue_raw == "" ~ "Unknown",
      tissue_lower == "rectum" ~ "Rectum",
      tissue_lower == "terminal_ileum" ~ "Terminal ileum",
      tissue_lower == "terminal ileum" ~ "Terminal ileum",
      tissue_lower == "ascending_colon" ~ "Ascending colon",
      tissue_lower == "ascending colon" ~ "Ascending colon",
      tissue_lower == "descending_colon" ~ "Descending colon",
      tissue_lower == "descending colon" ~ "Descending colon",
      tissue_lower == "sigmoid" ~ "Sigmoid colon",
      tissue_lower == "sigma" ~ "Sigmoid colon",
      tissue_lower == "transverse colon" ~ "Transverse colon",
      tissue_lower == "transverse_colon" ~ "Transverse colon",
      tissue_lower == "inflamed colon region" ~ "Inflamed colon region",
      tissue_lower == "healthy control colon" ~ "Healthy control colon",
      TRUE ~ tissue_raw
    )
  )

# 样本层面数据：用于 sex / age / disease 覆盖
sample_df <- meta_cov %>%
  dplyr::distinct(sample, disease, sex, age_num)

# specimen 层面数据：用于 tissue 覆盖
# 如果一个 sample 对应多个 tissue，这里会保留 sample-tissue 组合
specimen_df <- meta_cov %>%
  dplyr::distinct(sample, disease, sex, age_num, tissue_clean)

# =========================
# 2. 元数据完整性
# =========================

missing_df <- data.frame(
  variable = c("Disease", "Sex", "Age", "Tissue"),
  total = c(
    nrow(sample_df),
    nrow(sample_df),
    nrow(sample_df),
    nrow(specimen_df)
  ),
  missing = c(
    sum(is.na(sample_df$disease)),
    sum(is.na(sample_df$sex)),
    sum(is.na(sample_df$age_num)),
    sum(is.na(specimen_df$tissue_clean) | specimen_df$tissue_clean == "Unknown")
  )
) %>%
  dplyr::mutate(
    complete = total - missing,
    completeness = complete / total * 100,
    variable = factor(variable, levels = c("Disease", "Sex", "Age", "Tissue"))
  )

p_missing <- ggplot(missing_df, aes(x = variable, y = completeness)) +
  geom_col(width = 0.65, fill = "#6BAED6", color = "black", linewidth = 0.3) +
  geom_text(
    aes(label = paste0(round(completeness, 1), "%")),
    vjust = -0.4,
    size = 3.5
  ) +
  scale_y_continuous(
    limits = c(0, 105),
    breaks = seq(0, 100, 25),
    expand = c(0, 0)
  ) +
  labs(
    x = NULL,
    y = "Completeness (%)",
    title = "Metadata completeness"
  ) +
  theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    axis.text.x = element_text(color = "black"),
    axis.text.y = element_text(color = "black"),
    axis.title.y = element_text(face = "bold")
  )

# =========================
# 3. 性别覆盖
# =========================

sex_df <- sample_df %>%
  dplyr::filter(!is.na(disease), !is.na(sex)) %>%
  dplyr::count(disease, sex, name = "n")

p_sex <- ggplot(sex_df, aes(x = disease, y = n, fill = sex)) +
  geom_col(
    position = position_dodge(width = 0.7),
    width = 0.6,
    color = "black",
    linewidth = 0.3
  ) +
  geom_text(
    aes(label = n),
    position = position_dodge(width = 0.7),
    vjust = -0.4,
    size = 3.5
  ) +
  scale_fill_manual(values = c("Male" = "#96cada", "Female" = "#fdb26f")) +
  labs(
    x = NULL,
    y = "Number of samples",
    fill = "Sex",
    title = "Sex coverage"
  ) +
  theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    axis.text.x = element_text(color = "black", face = "bold"),
    axis.text.y = element_text(color = "black"),
    axis.title.y = element_text(face = "bold"),
    legend.title = element_text(face = "bold")
  )

# =========================
# 4. 年龄覆盖
# =========================

p_age <- ggplot(
  sample_df %>% dplyr::filter(!is.na(disease), !is.na(age_num)),
  aes(x = disease, y = age_num, fill = disease)
) +
  geom_violin(
    width = 0.75,
    alpha = 0.55,
    color = "black",
    linewidth = 0.3,
    trim = FALSE
  ) +
  geom_boxplot(
    width = 0.18,
    outlier.shape = NA,
    color = "black",
    linewidth = 0.35,
    alpha = 0.8
  ) +
  geom_jitter(
    width = 0.08,
    size = 1.8,
    alpha = 0.8
  ) +
  scale_fill_manual(values = c("HC" = "#96cada", "UC" = "#fdb26f")) +
  labs(
    x = NULL,
    y = "Age",
    title = "Age coverage"
  ) +
  theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    axis.text.x = element_text(color = "black", face = "bold"),
    axis.text.y = element_text(color = "black"),
    axis.title.y = element_text(face = "bold"),
    legend.position = "none"
  )

# =========================
# 5. Tissue 组成覆盖，按 disease 归一化
# =========================

tissue_df <- specimen_df %>%
  dplyr::filter(!is.na(disease), !is.na(tissue_clean), tissue_clean != "Unknown") %>%
  dplyr::count(disease, tissue_clean, name = "n") %>%
  dplyr::group_by(disease) %>%
  dplyr::mutate(prop = n / sum(n) * 100) %>%
  dplyr::ungroup()

tissue_order <- tissue_df %>%
  dplyr::group_by(tissue_clean) %>%
  dplyr::summarise(total_n = sum(n), .groups = "drop") %>%
  dplyr::arrange(dplyr::desc(total_n)) %>%
  dplyr::pull(tissue_clean)

tissue_df$tissue_clean <- factor(tissue_df$tissue_clean, levels = tissue_order)

tissue_colors <- setNames(
  hcl.colors(length(tissue_order), palette = "Dynamic"),
  tissue_order
)

p_tissue_stack <- ggplot(tissue_df, aes(x = disease, y = prop, fill = tissue_clean)) +
  geom_col(
    width = 0.6,
    color = "black",
    linewidth = 0.25
  ) +
  scale_fill_manual(values = tissue_colors) +
  scale_y_continuous(
    limits = c(0, 100),
    breaks = seq(0, 100, 25),
    expand = c(0, 0)
  ) +
  labs(
    x = NULL,
    y = "Tissue proportion (%)",
    fill = "Tissue",
    title = "Tissue coverage by disease group"
  ) +
  theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    axis.text.x = element_text(color = "black", face = "bold"),
    axis.text.y = element_text(color = "black"),
    axis.title.y = element_text(face = "bold"),
    legend.title = element_text(face = "bold"),
    legend.text = element_text(size = 8)
  )

# =========================
# 6. Disease × Tissue 覆盖热图
# =========================

heatmap_df <- specimen_df %>%
  dplyr::filter(!is.na(disease), !is.na(tissue_clean), tissue_clean != "Unknown") %>%
  dplyr::count(tissue_clean, disease, name = "n")

heatmap_df$tissue_clean <- factor(heatmap_df$tissue_clean, levels = rev(tissue_order))

p_tissue_heatmap <- ggplot(heatmap_df, aes(x = disease, y = tissue_clean, fill = n)) +
  geom_tile(color = "white", linewidth = 0.6) +
  geom_text(aes(label = n), size = 3.5, color = "black") +
  scale_fill_gradient(low = "#F7FBFF", high = "#2171B5") +
  labs(
    x = NULL,
    y = NULL,
    fill = "N",
    title = "Disease × tissue coverage"
  ) +
  theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    axis.text.x = element_text(color = "black", face = "bold"),
    axis.text.y = element_text(color = "black"),
    legend.title = element_text(face = "bold")
  )

# =========================
# 7. 合并展示
# =========================

coverage_plot <- (p_missing | p_sex | p_age) /
  (p_tissue_stack | p_tissue_heatmap) +
  plot_annotation(
    title = "Cohort metadata coverage",
    subtitle = "Coverage was summarized at sample or sample-tissue specimen level rather than cell level",
    theme = theme(
      plot.title = element_text(face = "bold", size = 17, hjust = 0.5),
      plot.subtitle = element_text(size = 11, hjust = 0.5, color = "grey35")
    )
  )

coverage_plot
