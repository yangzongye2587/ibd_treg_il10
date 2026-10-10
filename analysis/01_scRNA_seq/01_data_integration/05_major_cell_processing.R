cellRatioPlot(object = major_annotation_ibd_subtype_tregback,
              sample.name = "cell_subtype",
              celltype.name = "disease",fill.col = c('#96cada','#fdb26f'))+RotatedAxis()
# 运行原函数并赋值给p
p <- cellRatioPlot(object = major_annotation_ibd_subtype_tregback,
                   sample.name = "cell_subtype",
                   celltype.name = "disease",
                   fill.col = c('#96cada','#fdb26f'))

# 提取p中的绘图数据
df <- p$data

# 运行 View(df) 看一下数据结构，通常会包含几列：
# 一列是样本/亚型（比如 cell_subtype）
# 一列是分组（比如 disease）  
# 一列是比例值（比如 ratio, proportion 或 percentage）
View(df)
library(ggplot2)

# 如果是堆叠柱状图（通常细胞比例图是这种）：
ggplot(df, aes(x = cell_subtype, y = rel_num, fill = disease)) +
  # 将 color 从 "white" 改为 "black"，即可添加黑色描边
  geom_col(position = "stack", color = "black", width = 0.7) +
  scale_fill_manual(values = c('#96cada','#fdb26f')) +
  # 添加数字标签：position_stack(vjust=0.5)可以让数字显示在色块正中间
  geom_text(aes(label = ifelse(rel_num > 0.05, round(rel_num, 2), "")), 
            position = position_stack(vjust = 0.5), 
            size = 3, color = "black") +
  theme_classic() +
  labs(x = NULL, y = "Cell Ratio") + 
  RotatedAxis() # 注意：RotatedAxis()只需要写一次即可




cell_colors <- c(
  "#96cada",   # Light blue
  "#4eb144",   # Light green
  "#a7d492",          # Light gray
  "#347bae",           # Yello
  "#fb8484",      # Khaki
  "#fdb26f",           # Thistle
  "#edcc55",           # Light gray (similar to vγT)
  "#ff6060",        # Purple
  "#c3a9ca",        # Tomato
  "#8f76b2",           # Hot pink
  "#4c999f",     # Yellow (same as NK)
  '#f88517'
)
table(major_annotation_ibd_subtype_tregback$characteristics_ch1.14)
table(major_annotation_ibd_subtype_tregback$disease)
table(major_annotation_ibd_subtype_tregback$orig.ident)
library(Seurat)
library(ggplot2)
library(dplyr)
library(tidyr)

# ============================================
# 1. 提取每个病人的炎症评分（去除没有炎症评分的病人）
# ============================================

meta <- major_annotation_ibd_subtype_tregback@meta.data

# 提取炎症评分数值
meta$inflammation_score <- as.numeric(
  gsub("inflammation score: ", "", meta$characteristics_ch1.14)
)

# 按病人聚合，每个病人取一个炎症评分（同一病人评分相同）
patient_scores <- meta %>%
  filter(!is.na(inflammation_score)) %>%       # 去掉没有炎症评分的
  group_by(orig.ident, disease) %>%
  summarise(inflammation_score = first(inflammation_score), .groups = "drop")

cat("HC病人数:", sum(patient_scores$disease == "HC"), "\n")
cat("UC病人数:", sum(patient_scores$disease == "UC"), "\n")

# ============================================
# 2. 方法一：ggplot2 density 曲线（推荐）
# ============================================

p1 <- ggplot(patient_scores, aes(x = inflammation_score, fill = disease, color = disease)) +
  geom_density(alpha = 0.4, linewidth = 1.2) +
  scale_fill_manual(values = c("HC" = "#4DBBD5", "UC" = "#E64B35")) +
  scale_color_manual(values = c("HC" = "#4DBBD5", "UC" = "#E64B35")) +
  labs(
    title = "Inflammation Score Density by Disease",
    x = "Inflammation Score",
    y = "Density",
    fill = "Disease",
    color = "Disease"
  ) +
  theme_classic(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    legend.position = "top"
  )

ggsave("inflammation_density_curve.pdf", p1, width = 8, height = 6)
ggsave("inflammation_density_curve.png", p1, width = 8, height = 6, dpi = 300)
print(p1)

# ============================================
# 3. 方法二：带 rug plot 的增强版
# ============================================

p2 <- ggplot(patient_scores, aes(x = inflammation_score, fill = disease, color = disease)) +
  geom_density(alpha = 0.35, linewidth = 1.2) +
  geom_rug(aes(color = disease), alpha = 0.6, length = unit(0.08, "npc")) +
  scale_fill_manual(values = c("HC" = "#4DBBD5", "UC" = "#E64B35")) +
  scale_color_manual(values = c("HC" = "#4DBBD5", "UC" = "#E64B35")) +
  labs(
    title = "Inflammation Score Distribution: HC vs UC",
    subtitle = paste0("HC n=", sum(patient_scores$disease == "HC"),
                      " | UC n=", sum(patient_scores$disease == "UC")),
    x = "Inflammation Score",
    y = "Density"
  ) +
  theme_classic(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5, size = 11, color = "grey40"),
    legend.position = "top"
  )

ggsave("inflammation_density_rug.pdf", p2, width = 8, height = 6)

# ============================================
# 4. 方法三：分面显示（各自Y轴独立，更清晰）
# ============================================

p3 <- ggplot(patient_scores, aes(x = inflammation_score, fill = disease)) +
  geom_density(alpha = 0.6, color = "black", linewidth = 0.8) +
  facet_wrap(~ disease, scales = "free_y", nrow = 2) +
  scale_fill_manual(values = c("HC" = "#4DBBD5", "UC" = "#e9ad78")) +
  labs(
    title = "Inflammation Score by Disease (Faceted)",
    x = "Inflammation Score",
    y = "Density"
  ) +
  theme_classic(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    strip.text = element_text(face = "bold", size = 13)
  )
p3 
ggsave("inflammation_facet_density.pdf", p3, width = 8, height = 8)

# ============================================
# 5. 方法四：KDE平滑拟合 + 均值/中位数标注
# ============================================

stats_summary <- patient_scores %>%
  group_by(disease) %>%
  summarise(
    median_score = median(inflammation_score),
    mean_score   = mean(inflammation_score),
    n            = n()
  )

p4 <- ggplot(patient_scores, aes(x = inflammation_score, fill = disease, color = disease)) +
  geom_density(alpha = 0.35, linewidth = 1.2) +
  geom_vline(data = stats_summary, aes(xintercept = median_score, color = disease),
             linetype = "dashed", linewidth = 0.9) +
  annotate("text",
           x = stats_summary$median_score,
           y = c(0.05, 0.08),
           label = paste0("median=", round(stats_summary$median_score, 2)),
           color = c("#4DBBD5", "#e9ad78"),
           size = 4, fontface = "bold") +
  scale_fill_manual(values = c("HC" = "#4DBBD5", "UC" = "#e9ad78")) +
  scale_color_manual(values = c("HC" = "#4DBBD5", "UC" = "#e9ad78")) +
  labs(
    title = "Inflammation Score: HC vs UC (with Median)",
    x = "Inflammation Score",
    y = "Density"
  ) +
  theme_classic(base_size = 14) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"),
        legend.position = "top")

ggsave("inflammation_density_median.pdf", p4, width = 8, height = 6)

# ============================================
# 6. 打印统计摘要
# ============================================
cat("\n========== Inflammation Score Summary ==========\n")
print(stats_summary)
seurat_obj <- major_annotation_ibd_subtype_tregback
cell_order <- c('T/NK','plasma','B','endothelial_cells','Stromal',
                'Epithelial','MYLOID','pericyte','glia','other','Mast','neuron')

lpx_col <- setNames(
  c('#b47862','#F8862A','#FFBD6F','#727482','#ce2827',
    '#d4bcd7','#baa122','#d7e09d','#89982a','#359439','#6c95a3','#f19baa'),
  cell_order
)
# 年龄分布小提琴图
ggplot(major_annotation_ibd_subtype_tregback@meta.data, 
       aes(x = disease, y = age, fill = disease)) +
  geom_violin(alpha = 0.6, linewidth = 0.5) +
  geom_boxplot(width = 0.1, fill = "white", color = "black", 
               outlier.shape = NA, linewidth = 0.5) +
  scale_fill_manual(values = c("#6B9E8A", "#CB5B56")) +
  theme_classic(base_size = 8) +
  theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
        legend.position = "none") +
  labs(x = "Disease", y = "Age (years)")
# 性别分布堆叠柱状图
sex_dist <- major_annotation_ibd_subtype_tregback %>%
  group_by(disease, sex1) %>%
  summarise(count = n()) %>%
  mutate(prop = count / sum(count))

# 性别分布堆叠柱状图
sex_dist <- major_annotation_ibd_subtype_tregback@meta.data %>%
  group_by(disease, sex1) %>%
  summarise(count = n()) %>%
  mutate(prop = count / sum(count))

ggplot(sex_dist, aes(x = disease, y = prop, fill = sex1)) +
  geom_bar(stat = "identity", position = "fill", width = 0.7) +
  scale_fill_manual(values = c("#E8A0BF", "#8ABEDB")) +
  theme_classic(base_size = 8) +
  theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
        legend.position = "right",
        legend.title = element_blank()) +
  labs(x = "Disease", y = "Proportion") +
  scale_y_continuous(labels = scales::percent)
features <- c(
    "CD3E", 'MME', "IGHA2", "ACTA2",'SYT1', "LYZ",
  "KIT", "SOX10", "EPCAM", "PECAM1","CD19"
)

if (!"MKI67" %in% rownames(obj) && "Mki67" %in% rownames(obj)) {
  features[features == "MKI67"] <- "Mki67"
}

features_use <- features[features %in% rownames(obj)]

p <- DotPlot(
  object = obj,
  features = features_use,
  group.by = "cell_subtype"
)

dot_data <- p$data


p_dot_anno <- ggplot(dot_data, aes(x = features.plot, y = id)) +
  geom_point(
    aes(size = pct.exp, fill = avg.exp.scaled),
    shape = 21,
    color = "grey25",
    stroke = 0.25
  ) +
  scale_fill_gradient2(
    low = "#0072B2",
    mid = "white",
    high = "#D55E00",
    midpoint = 0,
    name = "Average\nexpression"
  ) +
  scale_size(
    range = c(0.2, 5),
    name = "Percent\nexpressed"
  ) +
  scale_x_discrete(position = "top") +
  labs(x = NULL, y = NULL) +
  theme_classic(base_size = 11) +
  theme(
    axis.text.x.top = element_text(
      angle = 45,
      hjust = 0,
      vjust = 0,
      color = "black",
      face = "bold",
      size = 9
    ),
    axis.text.y = element_text(
      color = "black",
      face = "bold",
      size = 9
    ),
    axis.ticks = element_blank(),
    axis.line = element_blank(),
    panel.grid.major = element_line(
      color = "grey90",
      linewidth = 0.25
    ),
    legend.position = "right",
    legend.title = element_text(face = "bold", size = 9),
    legend.text = element_text(size = 8)
  )

p_dot_anno




p_marker_heat <- ggplot(
  dot_data,
  aes(x = features.plot, y = id, fill = avg.exp.scaled)
) +
  geom_tile(color = "white", linewidth = 0.4) +
  scale_fill_gradient2(
    low = "#0072B2",
    mid = "white",
    high = "#D55E00",
    midpoint = 0,
    name = "Scaled\nexpression"
  ) +
  scale_x_discrete(position = "top") +
  labs(x = NULL, y = NULL) +
  theme_classic(base_size = 11) +
  theme(
    axis.text.x.top = element_text(
      angle = 45,
      hjust = 0,
      vjust = 0,
      color = "black",
      face = "bold",
      size = 9
    ),
    axis.text.y = element_text(
      color = "black",
      face = "bold",
      size = 9
    ),
    axis.ticks = element_blank(),
    axis.line = element_blank(),
    legend.position = "right"
  )

p_marker_heat

