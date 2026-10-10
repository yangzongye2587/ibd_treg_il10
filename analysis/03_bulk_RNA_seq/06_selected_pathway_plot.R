# ========== 开头：加载包、设置路径 ==========
setwd("data\\bulkil10\\过程文件\\JY")
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggnewscale)
library(stringr)
library(svglite)

dt1 <- read.csv("DEGs1.csv")
dt1 <- dt1 %>%
  mutate(
    group = str_squish(group),
    group = recode(
      group,
      "Type" = "IFN-gamma response"
    )
  )

# ============================================================
# 2. 清洗列名
# ============================================================
names(dt1) <- names(dt1) %>%
  str_replace_all("[\r\n]+", "_") %>%
  str_replace_all("\\s+", "_") %>%
  str_replace_all("\\.+", "_") %>%
  str_replace_all("_+", "_") %>%
  str_remove("^_") %>%
  str_remove("_$")

print(names(dt1))

# ============================================================
# 3. 检查所需列
# ============================================================
required_columns <- c(
  "gene_name",
  "group",
  "IN_vs_IBD_log2FoldChange",
  "IN_vs_IBD_pvalue",
  "TN_vs_IBD_log2FoldChange",
  "TN_vs_IBD_pvalue",
  "TN_vs_IN_log2FoldChange",
  "TN_vs_IN_pvalue"
)

missing_columns <- setdiff(required_columns, names(dt1))
if (length(missing_columns) > 0) {
  stop("以下列在数据中不存在：\n", paste(missing_columns, collapse = "\n"))
}

# ============================================================
# 4. 设置分组顺序、比较顺序
# ============================================================
group_order <- c(
  "Neutrophils",
  "Type I Interferon",
  "Th17",
  "Tissue injury"
)

compare_order <- c(
  "IN vs IBD",
  "TN vs IBD",
  "TN vs IN"
)

unknown_groups <- setdiff(unique(na.omit(dt1$group)), group_order)
if (length(unknown_groups) > 0) {
  warning("以下分组未包含在 group_order 中：", paste(unknown_groups, collapse = ", "))
}

gene_order <- dt1 %>%
  mutate(group = factor(group, levels = group_order)) %>%
  arrange(group, gene_name) %>%
  pull(gene_name) %>%
  unique()

# ============================================================
# 5. 宽表转长表
# ============================================================
dt_long <- dt1 %>%
  select(all_of(required_columns)) %>%
  pivot_longer(
    cols = -c(gene_name, group),
    names_to = c("comparison_code", "metric"),
    names_pattern = "(IN_vs_IBD|TN_vs_IBD|TN_vs_IN)_(log2FoldChange|pvalue)",
    values_to = "value"
  ) %>%
  pivot_wider(names_from = metric, values_from = value) %>%
  rename(log2FC = log2FoldChange, pval = pvalue) %>%
  mutate(
    compare = recode(
      comparison_code,
      "IN_vs_IBD" = "IN vs IBD",
      "TN_vs_IBD" = "TN vs IBD",
      "TN_vs_IN"  = "TN vs IN"
    ),
    y_position = case_when(
      compare == "IN vs IBD" ~ 3,
      compare == "TN vs IBD" ~ 2,
      compare == "TN vs IN"  ~ 1
    ),
    p_class = case_when(
      is.na(pval)  ~ "NS",
      pval < 0.01  ~ "p < 0.01",
      pval < 0.05  ~ "p < 0.05",
      TRUE         ~ "NS"
    ),
    p_class = factor(p_class, levels = c("NS", "p < 0.05", "p < 0.01")),
    gene_name = factor(gene_name, levels = gene_order),
    group = factor(group, levels = group_order)
  )

# ============================================================
# 6. 基因分组颜色
# ============================================================
group_colors <- c(
  "Neutrophils"       = "#D73027",
  "Type I Interferon" = "#4575B4",
  "Th17"              = "#FDAE61",
  "Tissue injury"     = "#984EA3"
)

# ============================================================
# 7. 绘图函数 (竖版：基因在Y轴，比较组在X轴) 固定log2FC -4 ~ 4
# ============================================================
gene_dot_plot_vertical <- function(data, title_text, custom_gene_order) {
  
  group_annotation <- data %>% distinct(gene_name, group)
  
  ggplot() +
    
    # 交替背景 (改为垂直方向交替)
    annotate("rect", ymin = -Inf, ymax = Inf, xmin = 2.5, xmax = 3.5,
             fill = "#F5F5F5", alpha = 0.8) +
    annotate("rect", ymin = -Inf, ymax = Inf, xmin = 0.5, xmax = 1.5,
             fill = "#F5F5F5", alpha = 0.8) +
    
    # 竖向分隔线 (分隔不同的比较组)
    geom_vline(xintercept = c(0.5, 1.5, 2.5, 3.5),
               color = "#D0D0D0", linewidth = 0.4) +
    
    # 横向网格线 (分隔不同的基因)
    geom_hline(yintercept = seq_along(custom_gene_order) + 0.5,
               color = "#E5E5E5", linewidth = 0.25) +
    
    # 左侧基因分组色条
    geom_tile(
      data = group_annotation,
      aes(x = 0.18, y = gene_name, fill = group),
      width = 0.18, height = 0.94, color = NA
    ) +
    scale_fill_manual(
      name = "Gene group",
      values = group_colors,
      drop = FALSE,
      na.value = "#BDBDBD"
    ) +
    guides(
      fill = guide_legend(
        order = 2, nrow = 2, byrow = TRUE,
        override.aes = list(linewidth = 0, alpha = 1)
      )
    ) +
    
    # 第二套 fill 色标
    ggnewscale::new_scale_fill() +
    
    # 气泡点 (x和y对调)
    geom_point(
      data = data,
      aes(x = y_position, y = gene_name, fill = log2FC,
          size = p_class, shape = p_class),
      color = "#404040", stroke = 0.5, alpha = 0.95, na.rm = TRUE
    ) +
    
    # ========== 固定尺度 -4 ~ 4 ==========
  scale_fill_gradient2(
    name = expression(log[2] * "(Fold Change)"),
    low = "#2166AC", mid = "#F7F7F7", high = "#B2182B",
    midpoint = 0,
    limits = c(-4, 4),
    oob = scales::squish,
    na.value = "#D9D9D9",
    breaks = c(-4,-2,0,2,4)
  ) +
    
    # P 值大小
    scale_size_manual(
      name = "P-value",
      values = c("NS" = 1.8, "p < 0.05" = 3, "p < 0.01" = 3.8),
      drop = FALSE
    ) +
    
    # P 值形状
    scale_shape_manual(
      name = "P-value",
      values = c("NS" = 24, "p < 0.05" = 21, "p < 0.01" = 21),
      drop = FALSE
    ) +
    
    # 纵坐标 (基因名，反转以使第一个基因在顶部)
    scale_y_discrete(
      limits = rev(custom_gene_order),
      drop = FALSE,
      expand = expansion(add = 0.55)
    ) +
    
    # 横坐标 (比较组)
    scale_x_continuous(
      breaks = c(3, 2, 1),
      labels = compare_order,
      limits = c(0.02, 3.5),
      expand = expansion(mult = c(0, 0))
    ) +
    
    # 标题
    labs(
      title = title_text,
      subtitle = paste0(
        "Color indicates log2 fold change; ",
        "shape and size indicate statistical significance"
      ),
      x = NULL, y = NULL,
      caption = paste0(
        "Triangles represent non‑significant results ",
        "(P ≥ 0.05 or missing P‑values)"
      )
    ) +
    
    # 图例
    guides(
      fill = guide_colorbar(
        order = 1, title.position = "top", title.hjust = 0.5,
        barwidth = grid::unit(4.3, "cm"),
        barheight = grid::unit(0.42, "cm")
      ),
      size = guide_legend(
        order = 2, title.position = "top",
        override.aes = list(
          shape = c(24, 21, 21),
          size = c(1.8, 3, 3.8),
          fill = "white",
          color = "#404040",
          alpha = 1, stroke = 0.5
        )
      ),
      shape = guide_legend(order = 2, title.position = "top")
    ) +
    
    coord_cartesian(clip = "off") +
    
    # 主题
    theme_classic(base_size = 11) +
    theme(
      plot.title = element_text(
        hjust = 0.5, face = "bold", size = 15,
        color = "#222222", margin = margin(b = 6)
      ),
      plot.subtitle = element_text(
        hjust = 0.5, size = 10,
        color = "#666666", margin = margin(b = 16)
      ),
      plot.caption = element_text(
        hjust = 1, size = 8.5,
        color = "#777777", margin = margin(t = 10)
      ),
      # X轴变为比较组，水平显示
      axis.text.x = element_text(
        size = 10, face = "bold",
        color = "#333333", margin = margin(t = 8)
      ),
      # Y轴变为基因名，斜体显示
      axis.text.y = element_text(
        size = 8.5, face = "italic", color = "#333333",
        margin = margin(r = 8)
      ),
      axis.ticks.x = element_blank(),
      axis.ticks.y = element_blank(),
      axis.line = element_blank(),
      panel.border = element_rect(
        color = "#555555", fill = NA, linewidth = 0.65
      ),
      legend.position = "bottom",
      legend.box = "vertical",
      legend.box.just = "center",
      legend.title = element_text(face = "bold", size = 9.5),
      legend.text = element_text(size = 8.5),
      legend.spacing.y = grid::unit(0.15, "cm"),
      plot.margin = margin(t = 16, r = 18, b = 12, l = 18)
    )
}

# ============================================================
# 8. 绘图
# ============================================================
p1 <- gene_dot_plot_vertical(
  data = dt_long,
  title_text = "Gene Expression Changes",
  custom_gene_order = gene_order
)

print(p1)
