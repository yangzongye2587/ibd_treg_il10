library(Seurat)
library(dplyr)
library(ggplot2)
library(stringr)
table(sce.all.filt$age)
sce.all.filt <- IBD_T_CELL_remove_dcs2
df <- sce.all.filt@meta.data %>%
  mutate(
    inflammation_score_num = as.numeric(str_extract(inflammation_score, "[0-9.]+$"))
  ) %>%
  filter(!is.na(age), !is.na(inflammation_score_num))

# 1) 线性拟合曲线
p_lm <- ggplot(df, aes(x = age, y = inflammation_score_num)) +
  geom_point(alpha = 0.25, size = 0.6) +
  geom_smooth(method = "lm", se = TRUE) +
  labs(x = "Age", y = "Inflammation score",
       title = "Inflammation score ~ Age (Linear fit)") +
  theme_classic()

print(p_lm)

# 2) LOESS 拟合（非参数，更像“平滑趋势线”）
p_loess <- ggplot(df, aes(x = age, y = inflammation_score_num)) +
  geom_point(alpha = 0.25, size = 0.6) +
  geom_smooth(method = "loess", se = TRUE, span = 0.75) +
  labs(x = "Age", y = "Inflammation score",
       title = "Inflammation score ~ Age (LOESS fit)") +
  theme_classic()

print(p_loess)

# 3) 二次项拟合（如果怀疑非线性但又想可解释）
p_poly2 <- ggplot(df, aes(x = age, y = inflammation_score_num,group = sex1)) +
  geom_point(alpha = 0.25, size = 0.6) +
  geom_smooth(method = "lm", formula = y ~ poly(x, 2, raw = TRUE), se = TRUE) +
  labs(x = "Age", y = "Inflammation score",
       title = "Inflammation score ~ Age (Quadratic fit)") +
  theme_classic()

print(p_poly2)

# 4) 相关性与模型统计（可选）
cat("Pearson cor:\n")
print(cor.test(df$age, df$inflammation_score_num, method = "pearson"))

cat("\nSpearman cor:\n")
print(cor.test(df$age, df$inflammation_score_num, method = "spearman"))

cat("\nLinear model summary:\n")
print(summary(lm(inflammation_score_num ~ age, data = df)))




library(ggplot2)
library(ggpmisc)

p_poly2 <- ggplot(df, aes(x = age, y = inflammation_score_num, color = sex1)) +
  geom_point(alpha = 0.25, size = 0.6) +
  geom_smooth(method = "lm", formula = y ~ poly(x, 2, raw = TRUE), se = TRUE) +
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..rr.label.., ..p.value.label.., sep = "~~~")),
    formula = y ~ poly(x, 2, raw = TRUE),
    parse = TRUE,
    label.x.npc = "left",   # 左侧
    label.y.npc = "top",    # 顶部
    size = 3
  ) +
  labs(x = "Age", y = "Inflammation score",
       title = "Inflammation score ~ Age (Quadratic fit by sex)") +
  theme_classic()

print(p_poly2)



neuron <- subset(ednrb_Stromal2.0,celltype_Stromal=='neuron'&time=='12d')
