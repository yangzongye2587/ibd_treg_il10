library(Seurat)
library(dplyr)
library(ggplot2)
library(patchwork)  # 拼图用

# ============================================================
# 0. 准备：将炎症评分字符串 → 数值（整个对象统一处理）
# ============================================================
IBD_T_TREG_SUB$inflammation_score_num <- as.numeric(
  sub("^inflammation score:\\s*",
      "",
      as.character(IBD_T_TREG_SUB$inflammation_score))
)

# ============================================================
# 1. 定义一个通用函数：给定 subset 后做 sample-level + cell-level 相关
# ============================================================

run_batf_cor <- function(seu, cell_type_label, cell_type_filter,
                         score_threshold = 2, batf_threshold = 0) {
  
  # ---- 取出目标亚型细胞 ----
  sub_cells <- subset(seu, subset = T_subtype == cell_type_filter)
  
  # ---- 提取 BATF + 炎症评分 + 样本 + 疾病 ----
  dat <- FetchData(sub_cells,
                   vars = c("BATF", "inflammation_score_num", "orig.ident", "disease"))
  colnames(dat)[colnames(dat) == "orig.ident"] <- "sample_id"
  
  # ============================
  # (A) Sample-level 相关（仅 UC）
  # ============================
  sample_uc <- dat %>%
    filter(disease == "UC" &
             !is.na(inflammation_score_num) &
             inflammation_score_num > score_threshold) %>%
    group_by(sample_id) %>%
    summarise(
      BATF_mean  = mean(BATF, na.rm = TRUE),
      inflam     = median(inflammation_score_num, na.rm = TRUE),
      n_cells    = n(),
      .groups    = "drop"
    ) %>%
    filter(!is.na(inflam))
  
  if (nrow(sample_uc) >= 3) {
    cor_sample <- cor.test(sample_uc$BATF_mean, sample_uc$inflam,
                           method = "spearman", exact = FALSE)
    rho_s <- round(unname(cor_sample$estimate), 3)
    p_s   <- signif(cor_sample$p.value, 3)
  } else {
    rho_s <- NA; p_s <- NA
  }
  
  p_sample <- ggplot(sample_uc, aes(x = inflam, y = BATF_mean)) +
    geom_point(size = 3, alpha = 0.8) +
    geom_smooth(method = "lm", se = TRUE, color = "steelblue", linewidth = 1) +
    annotate("text", x = Inf, y = Inf, hjust = 1.1, vjust = 1.5,
             label = sprintf("Spearman rho = %s\nP = %s\nSamples = %d",
                             rho_s, p_s, nrow(sample_uc)),
             size = 4.5) +
    theme_classic(base_size = 13) +
    labs(title    = paste0(cell_type_label, " — Sample-level (UC)"),
         subtitle = paste0("BATF mean vs Inflammation score (>", score_threshold, ")"),
         x = "Inflammation score", y = "Mean BATF expression")
  
  # ============================
  # (B) Cell-level 相关（仅 UC，BATF > 0 & score > threshold）
  # ============================
  cell_uc <- dat %>%
    filter(disease == "UC" &
             !is.na(inflammation_score_num) &
             inflammation_score_num > score_threshold &
             BATF > batf_threshold)
  
  if (nrow(cell_uc) >= 10) {
    cor_cell <- cor.test(cell_uc$BATF, cell_uc$inflammation_score_num,
                         method = "spearman", exact = FALSE)
    rho_c <- round(unname(cor_cell$estimate), 3)
    p_c   <- signif(cor_cell$p.value, 3)
  } else {
    rho_c <- NA; p_c <- NA
  }
  
  p_cell <- ggplot(cell_uc, aes(x = inflammation_score_num, y = BATF)) +
    geom_point(size = 0.5, alpha = 0.25) +
    geom_smooth(method = "lm", se = TRUE, color = "firebrick", linewidth = 1) +
    annotate("text", x = Inf, y = Inf, hjust = 1.1, vjust = 1.5,
             label = sprintf("Spearman rho = %s\nP = %s\nCells = %d",
                             rho_c, p_c, nrow(cell_uc)),
             size = 4.5) +
    theme_classic(base_size = 13) +
    labs(title    = paste0(cell_type_label, " — Cell-level (UC)"),
         subtitle = paste0("BATF > 0 & score >", score_threshold),
         x = "Inflammation score", y = "BATF expression")
  
  # 返回结果列表
  list(
    label       = cell_type_label,
    sample_df   = sample_uc,
    cor_sample  = cor_sample,
    cor_cell    = cor_cell,
    rho_sample  = rho_s,
    p_sample    = p_s,
    rho_cell    = rho_c,
    p_cell      = p_c,
    p_sample_plot = p_sample,
    p_cell_plot   = p_cell
  )
}

# ============================================================
# 2. 分别对 Treg / TH17 / TFH 运行
# ============================================================
res_treg <- run_batf_cor(IBD_T_TREG_SUB, "Treg",  "Treg")
res_th17 <- run_batf_cor(IBD_T_TREG_SUB, "TH17",  "TH17")
res_tfh  <- run_batf_cor(IBD_T_TREG_SUB, "TFH",   "TFH")

# 合并 sample-level 图
p_sample_all <- res_treg$p_sample_plot | res_th17$p_sample_plot | res_tfh$p_sample_plot
p_sample_all <- p_sample_all + plot_annotation(
  title = "Sample-level: BATF vs Inflammation score (UC, per CD4 subtype)",
  theme = theme(plot.title = element_text(size = 16, face = "bold", hjust = 0.5))
)

# 合并 cell-level 图
p_cell_all <- res_treg$p_cell_plot | res_th17$p_cell_plot | res_tfh$p_cell_plot
p_cell_all <- p_cell_all + plot_annotation(
  title = "Cell-level: BATF vs Inflammation score (UC, per CD4 subtype)",
  theme = theme(plot.title = element_text(size = 16, face = "bold", hjust = 0.5))
)

p_sample_all
p_cell_all

# ============================================================
# 3. CD4 总体（Treg + TH17 + TFH 合并）
# ============================================================
cd4_cells <- subset(IBD_T_TREG_SUB, subset = T_subtype %in% c("Treg", "TH17", "TFH"))

dat_cd4 <- FetchData(cd4_cells,
                     vars = c("BATF", "inflammation_score_num", "orig.ident", "disease"))
colnames(dat_cd4)[colnames(dat_cd4) == "orig.ident"] <- "sample_id"

# --- (A) Sample-level ---
sample_cd4_uc <- dat_cd4 %>%
  filter(disease == "UC" &
           !is.na(inflammation_score_num) &
           inflammation_score_num > 2) %>%
  group_by(sample_id) %>%
  summarise(
    BATF_mean = mean(BATF, na.rm = TRUE),
    inflam    = median(inflammation_score_num, na.rm = TRUE),
    n_cells   = n(),
    .groups   = "drop"
  ) %>%
  filter(!is.na(inflam))

cor_cd4_sample <- cor.test(sample_cd4_uc$BATF_mean, sample_cd4_uc$inflam,
                           method = "spearman", exact = FALSE)
rho_cd4_s <- round(unname(cor_cd4_sample$estimate), 3)
p_cd4_s   <- signif(cor_cd4_sample$p.value, 3)

p_cd4_sample <- ggplot(sample_cd4_uc, aes(x = inflam, y = BATF_mean)) +
  geom_point(size = 3, alpha = 0.8) +
  geom_smooth(method = "lm", se = TRUE, color = "darkgreen", linewidth = 1) +
  annotate("text", x = Inf, y = Inf, hjust = 1.1, vjust = 1.5,
           label = sprintf("Spearman rho = %s\nP = %s\nSamples = %d",
                           rho_cd4_s, p_cd4_s, nrow(sample_cd4_uc)),
           size = 5) +
  theme_classic(base_size = 14) +
  labs(title    = "CD4 total (Treg+TH17+TFH) — Sample-level (UC)",
       subtitle = "BATF mean vs Inflammation score (> 2)",
       x = "Inflammation score", y = "Mean BATF expression")

p_cd4_sample

# --- (B) Cell-level ---
cell_cd4_uc <- dat_cd4 %>%
  filter(disease == "UC" &
           !is.na(inflammation_score_num) &
           inflammation_score_num > 0.2 &
           BATF > 0)

cor_cd4_cell <- cor.test(cell_cd4_uc$BATF, cell_cd4_uc$inflammation_score_num,
                         method = "spearman", exact = FALSE)
rho_cd4_c <- round(unname(cor_cd4_cell$estimate), 3)
p_cd4_c   <- signif(cor_cd4_cell$p.value, 3)

p_cd4_cell <- ggplot(cell_cd4_uc, aes(x = inflammation_score_num, y = BATF)) +
  geom_point(size = 0.5, alpha = 0.2) +
  geom_smooth(method = "lm", se = TRUE, color = "purple", linewidth = 1) +
  annotate("text", x = Inf, y = Inf, hjust = 1.1, vjust = 1.5,
           label = sprintf("Spearman rho = %s\nP = %s\nCells = %d",
                           rho_cd4_c, p_cd4_c, nrow(cell_cd4_uc)),
           size = 5) +
  theme_classic(base_size = 14) +
  labs(title    = "CD4 total (Treg+TH17+TFH) — Cell-level (UC)",
       subtitle = "BATF > 0 & score > 0.2",
       x = "Inflammation score", y = "BATF expression")

p_cd4_cell

# ============================================================
# 4. 汇总表格
# ============================================================
summary_df <- data.frame(
  Group = c("Treg", "TH17", "TFH", "CD4_total"),
  `rho_sample` = c(res_treg$rho_sample, res_th17$rho_sample, res_tfh$rho_sample, rho_cd4_s),
  `p_sample`   = c(res_treg$p_sample,   res_th17$p_sample,   res_tfh$p_sample,   p_cd4_s),
  `rho_cell`   = c(res_treg$rho_cell,   res_th17$rho_cell,   res_tfh$rho_cell,   rho_cd4_c),
  `p_cell`     = c(res_treg$p_cell,     res_th17$p_cell,     res_tfh$p_cell,     p_cd4_c),
  check.names = FALSE
)
summary_df

# ============================================================
# 5. 把所有图拼在一起保存
# ============================================================
final_plot <- (p_sample_all / p_cell_all / (p_cd4_sample | p_cd4_cell)) +
  plot_annotation(
    title = "BATF vs Inflammation Score in CD4 T cell subtypes (UC)",
    theme = theme(plot.title = element_text(size = 18, face = "bold", hjust = 0.5))
  )

ggsave("BATF_inflammation_CD4_subtypes.pdf",
       final_plot, width = 18, height = 20, dpi = 300)
ggsave("BATF_inflammation_CD4_subtypes.png",
       final_plot, width = 18, height = 20, dpi = 300)
