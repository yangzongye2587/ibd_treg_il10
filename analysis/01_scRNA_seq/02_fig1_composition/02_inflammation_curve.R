library(Seurat)
library(survival)
library(survminer)
library(ggplot2)
library(dplyr)
library(tidyr)

# ============================================
# 1. 准备数据：每个病人一行
# ============================================
meta <- major_annotation_ibd_subtype_tregback@meta.data

meta$inflammation_score <- as.numeric(
  gsub("inflammation score: ", "", meta$characteristics_ch1.14)
)

patient_scores <- meta %>%
  filter(!is.na(inflammation_score)) %>%
  group_by(orig.ident, disease) %>%
  summarise(inflammation_score = first(inflammation_score), .groups = "drop")

# 构造生存分析所需变量
patient_scores$time   <- patient_scores$inflammation_score
patient_scores$status <- 1   # 全部为"事件"，无删失

cat("HC:", sum(patient_scores$disease == "HC"), "patients\n")
cat("UC:", sum(patient_scores$disease == "UC"), "patients\n")

# ============================================
# 2. 方法一：survminer 绘制经典 KM 曲线（推荐）
# ============================================
fit <- survfit(Surv(time, status) ~ disease, data = patient_scores)

# Log-rank 检验
logrank_p <- surv_pvalue(fit)$pval

p_km <- ggsurvplot(
  fit,
  data = patient_scores,
  
  # ---- 曲线样式 ----
  palette          = c("#4DBBD5", "#E64B35"),   # HC蓝, UC红
  risk.table       = TRUE,                      # 底部风险表
  risk.table.height = 0.25,
  risk.table.title = "Number at risk",
  
  # ---- 坐标轴 ----
  xlab = "Inflammation Score",
  ylab = "Proportion with score > x",
  xlim = c(0, max(patient_scores$time, na.rm = TRUE) + 0.5),
  
  # ---- 标题与P值 ----
  title    = "Inflammation Score: HC vs UC",
  pval     = TRUE,
  pval.method = TRUE,
  pval.size   = 5,
  
  # ---- 曲线细节 ----
  conf.int       = FALSE,
  surv.median.line = "hv",       # 中位线
  censor         = FALSE,        # 无删失点，不用画
  legend.title   = "Disease",
  legend.labs    = c("HC", "UC"),
  ggtheme = theme_classic(base_size = 14) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"))
)

pdf("inflammation_KM_curve.pdf", width = 9, height = 8)
print(p_km)
dev.off()

png("inflammation_KM_curve.png", width = 9, height = 8, units = "in", res = 300)
print(p_km)
dev.off()

# ============================================
# 3. 方法二：手动 ggplot2 绘制（更灵活定制）
# ============================================
# 计算每个 disease 组的 1-ECDF
ecdf_data <- patient_scores %>%
  arrange(disease, inflammation_score) %>%
  group_by(disease) %>%
  mutate(
    ecdf_val   = row_number() / n(),                # P(score ≤ x)
    surv_val   = 1 - ecdf_val,                      # P(score > x) ← 生存曲线
    n_at_risk  = n() - row_number() + 1
  ) %>%
  ungroup()

# 中位数炎症评分（曲线降到0.5的位置）
median_scores <- patient_scores %>%
  group_by(disease) %>%
  summarise(median_score = median(inflammation_score), .groups = "drop")

# KS 检验
ks_result <- ks.test(
  patient_scores$inflammation_score[patient_scores$disease == "HC"],
  patient_scores$inflammation_score[patient_scores$disease == "UC"]
)
ks_p_label <- paste0("KS test P = ", signif(ks_result$p.value, 3))

p_manual <- ggplot(ecdf_data, aes(x = inflammation_score, y = surv_val,
                                  color = disease, fill = disease)) +
  # 阶梯曲线
  geom_step(linewidth = 1.2, direction = "hv") +
  
  # 中位线
  geom_hline(yintercept = 0.5, linetype = "dashed", color = "grey50", linewidth = 0.5) +
  geom_vline(data = median_scores,
             aes(xintercept = median_score, color = disease),
             linetype = "dotted", linewidth = 0.8, show.legend = FALSE) +
  
  # 标注中位数
  geom_text(data = median_scores,
            aes(x = median_score, y = 0.57,
                label = paste0("median = ", round(median_score, 2)),
                color = disease),
            hjust = -0.05, size = 4, fontface = "bold", show.legend = FALSE) +
  
  # P值标注
  annotate("text", x = max(ecdf_data$inflammation_score) * 0.6, y = 0.15,
           label = ks_p_label, size = 5, fontface = "bold", color = "black") +
  
  # 颜色
  scale_color_manual(values = c("HC" = "#4DBBD5", "UC" = "#E64B35")) +
  scale_fill_manual(values  = c("HC" = "#4DBBD5", "UC" = "#E64B35")) +
  
  # 坐标轴
  scale_y_continuous(limits = c(0, 1.05), expand = c(0, 0)) +
  labs(
    title = "Inflammation Score: HC vs UC",
    subtitle = paste0("HC n = ", sum(patient_scores$disease == "HC"),
                      "  |  UC n = ", sum(patient_scores$disease == "UC")),
    x = "Inflammation Score",
    y = "Proportion with score > x",
    color = "Disease"
  ) +
  
  theme_classic(base_size = 14) +
  theme(
    plot.title    = element_text(hjust = 0.5, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5, size = 11, color = "grey40"),
    legend.position = "top"
  )

ggsave("inflammation_survival_style_curve.pdf", p_manual, width = 9, height = 7)
ggsave("inflammation_survival_style_curve.png", p_manual, width = 9, height = 7, dpi = 300)

# ============================================
# 4. 方法三：反向版本 — 累积发病曲线（ECDF 上升型）
# ============================================
p_ecdf <- ggplot(ecdf_data, aes(x = inflammation_score, y = ecdf_val,
                                color = disease, fill = disease)) +
  geom_step(linewidth = 1.2, direction = "hv") +
  geom_hline(yintercept = 0.5, linetype = "dashed", color = "grey50") +
  scale_color_manual(values = c("HC" = "#4DBBD5", "UC" = "#E64B35")) +
  scale_y_continuous(limits = c(0, 1.05), expand = c(0, 0)) +
  labs(
    title = "Cumulative Distribution of Inflammation Score",
    x = "Inflammation Score",
    y = "Proportion with score ≤ x",
    color = "Disease"
  ) +
  theme_classic(base_size = 14) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"),
        legend.position = "top")

ggsave("inflammation_ECDF_curve.pdf", p_ecdf, width = 9, height = 7)

# ============================================
# 5. 输出统计摘要
# ============================================
cat("\n========== Summary ==========\n")
cat("HC median inflammation score:",
    median(patient_scores$inflammation_score[patient_scores$disease == "HC"]), "\n")
cat("UC median inflammation score:",
    median(patient_scores$inflammation_score[patient_scores$disease == "UC"]), "\n")
cat("Log-rank P value:", signif(logrank_p, 4), "\n")
cat("KS test D =", signif(ks_result$statistic, 4),
    ", P =", signif(ks_result$p.value, 4), "\n")
