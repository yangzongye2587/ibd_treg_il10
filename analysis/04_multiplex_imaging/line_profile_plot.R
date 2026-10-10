# ========== 开头：仅加载包、设置路径 ==========
setwd("E:/Y/T")

library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)
library(ggnewscale)
library(stringr)
library(svglite)

dt1 <- read.csv("DEGs.csv")

# 通路分组顺序
group_order <- c(
  "Neutrophils",
  "Inflammatory factors",
  "Type I Interferon",
  "Type II Interferon",
  "Th17",
  "Tissue injury",
  "Homeostasis",
  "Treg"
)

gene_order <- dt1 %>%
  mutate(group = factor(group, levels = group_order)) %>%
  arrange(group, gene_name) %>%
  pull(gene_name) %>%
  unique()

#===== 宽表中分别对三组独立计算p_class，逻辑写三遍 =====
dt_part1 <- dt1 %>%
  select(gene_name, group, log2FC = IN_vs_IBD._log2FoldChange, pval = IN_vs_IBD._pvalue) %>%
  mutate(compare = "IN vs IBD") %>%
  mutate(
    p_class = case_when(
      is.na(pval) | pval > 0.05 ~ "NS",
      pval < 0.01 ~ "p<0.01",
      pval < 0.05 ~ "p<0.05"
    ),
    p_class = factor(p_class, levels = c("NS", "p<0.05", "p<0.01"))
  )

dt_part2 <- dt1 %>%
  select(gene_name, group, log2FC = TN_vs_IBD._log2FoldChange, pval = TN_vs_IBD._pvalue) %>%
  mutate(compare = "TN vs IBD") %>%
  mutate(
    p_class = case_when(
      is.na(pval) | pval > 0.05 ~ "NS",
      pval < 0.01 ~ "p<0.01",
      pval < 0.05 ~ "p<0.05"
    ),
    p_class = factor(p_class, levels = c("NS", "p<0.05", "p<0.01"))
  )

dt_part3 <- dt1 %>%
  select(gene_name, group, log2FC = TN_vs_IN_log2FoldChange, pval = TN_vs_IN_pvalue) %>%
  mutate(compare = "TN vs IN") %>%
  mutate(
    p_class = case_when(
      is.na(pval) | pval > 0.05 ~ "NS",
      pval < 0.01 ~ "p<0.01",
      pval < 0.05 ~ "p<0.05"
    ),
    p_class = factor(p_class, levels = c("NS", "p<0.05", "p<0.01"))
  )

# 拼接成长表
dt_long <- bind_rows(dt_part1, dt_part2, dt_part3) %>%
  mutate(
    gene_name = factor(gene_name, levels = gene_order),
    group = factor(group, levels = group_order)
  )

compare_order <- c("IN vs IBD", "TN vs IBD", "TN vs IN")

# ===================== 横向点图函数 =====================
gene_dot_plot_horizontal <- function(data, title_text) {
  
  fc_lim <- max(abs(data$log2FC), na.rm = TRUE)
  
  ggplot(data, aes(x = gene_name, y = compare)) +
    geom_rect(
      xmin = -Inf, xmax = Inf,
      ymin = 0.5, ymax = length(compare_order) + 0.5,
      fill = NA, color = "grey40", linewidth = 0.6,
      inherit.aes = FALSE
    ) +
    geom_hline(
      yintercept = seq_along(compare_order),
      color = "grey80", linewidth = 0.3, alpha = 0.7
    ) +
    geom_vline(
      xintercept = seq_along(gene_order),
      color = "grey80", linewidth = 0.3, alpha = 0.7
    ) +
    
    geom_tile(
      data = data %>% distinct(gene_name, group),
      aes(x = gene_name, y = 0.25, fill = group),
      width = 0.85, height = 0.22, alpha = 0.92,
      inherit.aes = FALSE
    ) +
    scale_fill_manual(
      name = "Gene group",
      values = setNames(
        c("#E64B35", "#4DBBD5", "#3C5488", "#F39B7F", "#8491B4", "#91D1C2", "#7E6148", "#B09C85"),
        group_order
      )
    ) +
    new_scale_fill() +
    
    geom_point(
      aes(fill = log2FC, size = p_class),
      shape = 21, stroke = 0.35, alpha = 0.88, color = "grey30"
    ) +
    
    scale_fill_gradient2(
      name = "log2(Fold Change)",
      low = "#2166ac", mid = "#f7f7f7", high = "#b2182b",
      midpoint = 0, limits = c(-fc_lim, fc_lim)
    ) +
    
    scale_size_manual(
      name = "P‑value",
      values = c("NS" = 1.4, "p<0.05" = 3.2, "p<0.01" = 5)
    ) +
    
    scale_x_discrete(limits = gene_order) +
    scale_y_discrete(limits = compare_order) +
    
    labs(x = NULL, y = NULL, title = title_text) +
    
    guides(
      fill = guide_colorbar(order = 1),
      size = guide_legend(order = 2)
    ) +
    
    theme_minimal(base_size = 10) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 12),
      axis.text.x = element_text(angle = 45, hjust = 1, size = 8.5),
      axis.text.y = element_text(size = 9, face = "bold"),
      panel.grid = element_blank(),
      plot.margin = margin(15,15,15,15)
    )
}

p_main_horizontal <- gene_dot_plot_horizontal(dt_long, "Gene expression changes across comparisons")
print(p_main_horizontal)

ggsave("gene_dotplot_horizontal.pdf", p_main_horizontal, width = 18, height = 9, dpi = 300)
ggsave("gene_dotplot_horizontal.svg", p_main_horizontal, width = 18, height = 9)