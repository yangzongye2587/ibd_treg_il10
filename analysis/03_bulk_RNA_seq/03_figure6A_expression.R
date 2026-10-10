# ========== 开头：加载包、设置路径 ==========
setwd("data\\bulkil10\\过程文件")
library(dplyr)
library(tidyr)
library(ggplot2)
library(stringr)

# ============================================================
# 1. 读取 FPKM 数据
# ============================================================
file_path <- "gene_sample_FPKM.xls"

dt_fpkm <- tryCatch(
  read.table(file_path, sep = "\t", header = TRUE, check.names = FALSE, quote = ""),
  error = function(e) {
    readxl::read_xls(file_path)
  }
)

names(dt_fpkm) <- names(dt_fpkm) %>%
  str_replace_all("[\r\n]+", "_") %>%
  str_replace_all("\\s+", "_") %>%
  str_replace_all("\\.+", "_") %>%
  str_replace_all("_+", "_") %>%
  str_remove("^_") %>%
  str_remove("_$")

# ============================================================
# 2. 提取目标基因 (Il10 和 Actb)
# ============================================================
hk_gene <- "Actb"
if (!hk_gene %in% dt_fpkm$gene_name) {
  stop("未在数据中找到内参基因 Actb，请检查 gene_name 列。")
}

target_genes <- c("Il10", hk_gene)
sample_cols <- c("IBD_1","IBD_2","IBD_3","TN_1","TN_2","TN_3","IN_1","IN_2","IN_3")

dt_target <- dt_fpkm %>%
  filter(gene_name %in% target_genes) %>%
  select(gene_name, all_of(sample_cols)) %>%
  pivot_longer(cols = all_of(sample_cols), names_to = "Sample", values_to = "FPKM") %>%
  mutate(
    Group = str_extract(Sample, "^[A-Z]+"),
    Group = factor(Group, levels = c("IBD", "TN", "IN")),
    gene_name = factor(gene_name, levels = c("Il10", hk_gene))
  )

# ============================================================
# 3. 计算均值 + 标准误
# ============================================================
dt_summary <- dt_target %>%
  group_by(Group, gene_name) %>%
  summarise(
    Mean = mean(FPKM, na.rm = TRUE),
    SE   = sd(FPKM, na.rm = TRUE) / sqrt(n()),
    .groups = "drop"
  )

# ============================================================
# 4. 配色与参数
# ============================================================
gene_colors <- c(
  "Il10" = "#D73027",
  "Actb" = "#4575B4"
)

dodge_width <- 0.8
bar_width <- 0.65

# ============================================================
# 5. 绘图
# ============================================================
p <- ggplot() +
  
  # ---- 交替背景 ----
annotate("rect", xmin = 0.5, xmax = 1.5, ymin = -Inf, ymax = Inf,
         fill = "#F5F5F5", alpha = 0.5) +
  annotate("rect", xmin = 2.5, xmax = 3.5, ymin = -Inf, ymax = Inf,
           fill = "#F5F5F5", alpha = 0.5) +
  
  # ---- 分隔线 ----
geom_vline(xintercept = c(1.5, 2.5), color = "#D0D0D0", linewidth = 0.3) +
  
  # ---- 柱子 (均值) ----
geom_col(
  data = dt_summary,
  aes(x = Group, y = Mean, fill = gene_name),
  position = position_dodge(width = dodge_width),
  width = bar_width,
  alpha = 0.85, color = NA
) +
  
  # ---- 误差棒 (标准误) ----
geom_errorbar(
  data = dt_summary,
  aes(x = Group, ymin = Mean - SE, ymax = Mean + SE, group = gene_name),
  position = position_dodge(width = dodge_width),
  width = 0.2, color = "#404040", linewidth = 0.5
) +
  
  # ---- 散点 (每个样本的真实值) ----
geom_point(
  data = dt_target,
  aes(x = Group, y = FPKM, fill = gene_name, shape = gene_name),
  position = position_dodge(width = dodge_width),
  color = "#333333", size = 2.5, 
  stroke = 0.6, alpha = 0.9
) +
  
  # ---- 填充色 (控制柱子和散点内部填充) ----
scale_fill_manual(
  name = "Gene",
  values = gene_colors,
  drop = FALSE,
  guide = guide_legend(
    override.aes = list(
      shape = c(21, 24),
      size = 4,
      color = "#333333",
      stroke = 0.6
    )
  )
) +
  
  # ---- 散点形状 ----
scale_shape_manual(
  name = "Gene",
  values = c("Il10" = 21, "Actb" = 24),
  drop = FALSE,
  guide = "none" # 隐藏 shape 的独立图例，合并到 fill 中
) +
  
  # ---- Y 轴 ----
scale_y_continuous(
  expand = expansion(mult = c(0, 0.08))
) +
  
  # ---- 标题 ----
labs(
  title = "Validation of Il10 Knockout",
  subtitle = "Il10 is ablated across all groups; housekeeping gene remains intact",
  x = NULL,
  y = "FPKM",
  caption = "Bars: mean ± SE; dots: individual samples (n=3 per group)"
) +
  
  # ---- 主题 ----
theme_classic(base_size = 12) +
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
    axis.text.x = element_text(
      size = 12, face = "bold", color = "#333333",
      margin = margin(t = 8)
    ),
    axis.text.y = element_text(
      size = 10, color = "#333333",
      margin = margin(r = 8)
    ),
    axis.title.y = element_text(
      size = 11, face = "bold", color = "#333333",
      margin = margin(r = 10)
    ),
    axis.ticks.x = element_blank(),
    axis.ticks.length.y = unit(0.15, "cm"),
    axis.line = element_line(color = "#555555", linewidth = 0.5),
    panel.border = element_rect(
      color = "#555555", fill = NA, linewidth = 0.65
    ),
    legend.position = "bottom",
    legend.title = element_text(face = "bold", size = 9.5),
    legend.text = element_text(size = 8.5, face = "italic"),
    legend.key.size = unit(0.5, "cm"),
    plot.margin = margin(t = 16, r = 16, b = 10, l = 14)
  )

print(p)
