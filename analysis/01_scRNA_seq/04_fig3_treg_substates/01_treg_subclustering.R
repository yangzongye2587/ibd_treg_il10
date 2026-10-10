setwd("data")  
AAA <- as.data.frame(rownames())
unique(TREG_SUBSET@meta.data$treg_sub)
TREG_SUBSET@meta.data$treg_sub[
  TREG_SUBSET@meta.data$treg_sub == "cyclin_TREG"
] <- "NAIVE_LIKE_TREG"
saveRDS(TREG_SUBSET,file = 'TREG_SUBSET.rds')
sce.all.filt <- TREG_SUBSET
sce.all.filt <- sce.all.filt <- sce.all.filt <- subset(IBD_T_CELL_remove_dcs2,T_subtype=='Treg')
sce.all.filt <- NormalizeData(sce.all.filt, 
                              normalization.method = "LogNormalize",
                              scale.factor = 1e4) 
sce.all.filt <- FindVariableFeatures(sce.all.filt)
sce.all.filt <- ScaleData(sce.all.filt)
sce.all.filt <- RunPCA(sce.all.filt, features = VariableFeatures(object = sce.all.filt))
library(harmony)
sce.all.filt <- RunHarmony(sce.all.filt, "orig.ident")

sce.all.filt <- RunUMAP(sce.all.filt,  dims = 1:15, 
                        reduction = "harmony")
DimPlot(sce.all.filt,reduction = "umap",label=T,group.by = 'orig.ident' ) 
sce.all.filt <- FindNeighbors(sce.all.filt, reduction = "harmony",
                              dims = 1:15) 
for (res in c(0.01, 0.05, 0.1, 0.2, 0.3, 0.5,0.8,1)) {
  sce.all.filt=FindClusters(sce.all.filt, #graph.name = "CCA_snn", 
                            resolution = res, algorithm = 1)
}
library(clustree)
p2_tree=clustree(sce.all.filt@meta.data, prefix = "RNA_snn_res.")
p2_tree
DimPlot(sce.all.filt,reduction = "umap",label=T,group.by = 'RNA_snn_res.0.2' ) 
TREG_SUBSET 
DotPlot(major_annotation_ibd_subtype_tregback,
        features =c('BATF','JAK1','JAK3','STAT3','STAT5','STAT1','IL2RA','TNFRSF9','TNFRSF4','TNFRSF1B','TNFRSF18','TNF','FOXP3','SELL','CCR7','FOS','JUN','DUSP1','MAF','IL10','IL22'),
        group.by = "final_subtype",cols = c('white','red'))+RotatedAxis()

cellRatioPlot(object = sce.all.filt,
              sample.name = "disease",
              celltype.name = "treg_sub")
sel.clust = "treg_sub"
sce.all.filt <- SetIdent(sce.all.filt, value = sel.clust)
necp0_markers<- FindAllMarkers(sce.all.filt, only.pos = TRUE, min.pct = 0.5, logfc.threshold = 0.25 )

celltype=data.frame(ClusterID=0:4,
                    celltype= 0:4) 
celltype[celltype$ClusterID %in% c(0),2]='NAIVE_LIKE_TREG'
celltype[celltype$ClusterID %in% c(1),2]='actived_TREG'
celltype[celltype$ClusterID %in% c(2),2]="actived_TREG"
celltype[celltype$ClusterID %in% c(3),2]='TNF_TREG'
celltype[celltype$ClusterID %in% c(4),2]='NAIVE_LIKE_TREG'
sce.all.filt@meta.data$treg_sub = "NA"
for(i in 1:nrow(celltype)){
  sce.all.filt@meta.data[which(sce.all.filt@meta.data$RNA_snn_res.0.2 == celltype$ClusterID[i]),'treg_sub'] <- celltype$celltype[i]}
sce.all.filt <- TREG_SUBSET
lpx_col <- c('#F8862A','#ce2827','#359439','#6c95a3','#b47862','#FFBD6F','#727482','#d4bcd7','#baa122',
             '#d7e09d','#89982a','#f19baa'
             ,'#6859a0','#a1d198','#678759')
DimPlot(sce.all.filt,reduction = "umap",label=T,group.by = 'treg_sub',cols = lpx_col ) 

cellRatioPlot(object = sce.all.filt,
              sample.name = "disease",
              celltype.name = "treg_sub",fill.col = lpx_col)
getwd()



saveRDS(sce.all.filt,file = 'TREG_SUBSET.RDS')
sce.all.filt <- TREG_SUBSET
library(Seurat)
library(dplyr)
library(tidyr)
library(purrr)
meta_df <- TREG_SUBSET@meta.data %>%
  mutate(
    cell = rownames(.),
    inflammation_score_num = 
      as.numeric(gsub("inflammation score: ", "", as.character(inflammation_score)))
    
  ) 
TREG_SUBSET@meta.data <- meta_df
head(TREG_SUBSET@meta.data)
library(dplyr)
library(tidyr)
library(ggplot2)

meta <- TREG_SUBSET@meta.data

# 只保留需要的列
df <- meta %>%
  dplyr::select(orig.ident, treg_sub, inflammation_score_num)

# 先看每个病人的 inflammation_score_num 是否唯一
check_score <- df %>%
  group_by(orig.ident) %>%
  summarise(n_score = n_distinct(inflammation_score_num), .groups = "drop")

table(check_score$n_score)
# 每个病人每个 treg_sub 的细胞数
sub_count <- df %>%
  group_by(orig.ident, treg_sub) %>%
  summarise(n_sub = n(), .groups = "drop")

# 每个病人的总 Treg 细胞数
patient_total <- df %>%
  group_by(orig.ident) %>%
  summarise(n_total = n(),
            inflammation_score_num = unique(inflammation_score_num)[1],
            .groups = "drop")

# 合并并计算比例
patient_prop_long <- sub_count %>%
  left_join(patient_total, by = "orig.ident") %>%
  mutate(prop = n_sub / n_total)

head(patient_prop_long)
cor_res <- patient_prop_long %>% filter(!inflammation_score_num == '') %>%
  group_by(treg_sub) %>%
  summarise(
    cor = suppressWarnings(cor(prop, inflammation_score_num, method = "spearman")),
    p = suppressWarnings(cor.test(prop, inflammation_score_num, method = "spearman")$p.value),
    .groups = "drop"
  ) %>%
  mutate(p_adj = p.adjust(p, method = "BH")) %>%
  arrange(p_adj)

cor_res
ggplot(patient_prop_long, aes(x = prop, y = inflammation_score_num)) +
  geom_point(size = 2, alpha = 0.8) +
  geom_smooth(method = "lm", se = TRUE) +
  facet_wrap(~treg_sub, scales = "free_x") +
  theme_bw() +
  labs(x = "Proportion of treg_sub in each patient",
       y = "Inflammation score",
       title = "Association between Treg subtype proportion and inflammation score")


library(dplyr)
library(tidyr)
library(ggplot2)

meta <- TREG_SUBSET@meta.data

df <- meta %>%
  dplyr::select(orig.ident, treg_sub, inflammation_score_num)

# 检查每个病人的 inflammation score 是否唯一
check_score <- df %>%
  group_by(orig.ident) %>%
  summarise(n_score = n_distinct(inflammation_score_num), .groups = "drop")

print(table(check_score$n_score))

# 病人-亚群计数
sub_count <- df %>%
  group_by(orig.ident, treg_sub) %>%
  summarise(n_sub = n(), .groups = "drop")

# 病人总数 + 炎症分数
patient_total <- df %>%
  group_by(orig.ident) %>%
  summarise(
    n_total = n(),
    inflammation_score_num = unique(inflammation_score_num)[1],
    .groups = "drop"
  )

# 比例表
patient_prop_long <- sub_count %>%
  left_join(patient_total, by = "orig.ident") %>%
  mutate(prop = n_sub / n_total)

# 相关分析
cor_res <- patient_prop_long %>%
  group_by(treg_sub) %>%
  summarise(
    cor = suppressWarnings(cor(prop, inflammation_score_num, method = "spearman")),
    p = suppressWarnings(cor.test(prop, inflammation_score_num, method = "spearman")$p.value),
    .groups = "drop"
  ) %>%
  mutate(p_adj = p.adjust(p, method = "BH")) %>%
  arrange(p_adj)

print(cor_res)

# 散点图
p1 <- ggplot(patient_prop_long, aes(x = prop, y = inflammation_score_num)) +
  geom_point(size = 2, alpha = 0.8) +
  geom_smooth(method = "lm", se = TRUE) +
  facet_wrap(~treg_sub, scales = "free_x") +
  theme_bw() +
  labs(x = "Proportion of treg_sub in each patient",
       y = "Inflammation score",
       title = "Treg subtype proportion vs inflammation score")

print(p1)



library(ggplot2)
library(ggpubr)

p1 <- ggplot(patient_prop_long, aes(x = prop, y = inflammation_score_num)) +
  geom_point(
    size = 2.5,
    alpha = 0.85,
    color = "#2C7FB8"
  ) +
  geom_smooth(
    method = "lm",
    se = TRUE,
    color = "#D95F0E",
    fill = "#FDD0A2",
    linewidth = 0.9
  ) +
  stat_cor(
    method = "spearman",
    label.x.npc = "left",
    label.y.npc = "top",
    aes(label = paste(..r.label.., ..p.label.., sep = "~`,`~")),
    size = 4
  ) +
  facet_wrap(~treg_sub, scales = "free_x") +
  labs(
    x = "Proportion of treg_sub in each patient",
    y = "Inflammation score",
    title = "Treg subtype proportion vs inflammation score"
  ) +
  theme_bw(base_size = 13) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
    axis.title = element_text(face = "bold"),
    axis.text = element_text(color = "black"),
    strip.text = element_text(face = "bold", size = 11),
    strip.background = element_rect(fill = "#F5F5F5", color = "black"),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(color = "grey90"),
    panel.border = element_rect(color = "black", linewidth = 0.8),
    plot.margin = margin(10, 15, 10, 15)
  )

print(p1)
# 1. 提取子集注释
sub_anno <- TREG_SUBSET@meta.data[, "treg_sub", drop = FALSE]
colnames(sub_anno) <- "treg_sub"
Tcell <- IBD_T_CELL_remove_dcs2
# 2. 回填到大群对象

Tcell <- AddMetaData(Tcell, metadata = sub_anno)

# 3. 可选：生成一个融合列
Tcell$T_subtype_refine <- Tcell$T_subtype
idx <- rownames(TREG_SUBSET@meta.data)
Tcell@meta.data[idx, "T_subtype_refine"] <- Tcell@meta.data[idx, "treg_sub"]

# 4. 检查
table(Tcell$T_subtype_refine, useNA = "ifany")
DimPlot(Tcell,group.by = 'T_subtype_refine')
head(Tcell)

library(dplyr)
library(tidyr)

library(dplyr)
library(tidyr)
library(stringr)
meta_df <- Tcell@meta.data %>%
  mutate(
    cell = rownames(.),
    inflammation_score_num = 
      as.numeric(gsub("inflammation score: ", "", as.character(inflammation_score)))
    
  ) 
Tcell@meta.data <- meta_df
meta <- Tcell@meta.data

# 提取病人级炎症信息
patient_info <- meta %>%
  group_by(orig.ident) %>%
  summarise(
    inflammation_score_num = unique(inflammation_score_num)[1],
    inflamed = unique(characteristics_ch1.7)[1],
    .groups = "drop"
  ) %>%
  mutate(
    inflamed = ifelse(str_detect(inflamed, "Inflamed"), 1, 0)
  )

head(patient_info)

target_sub <- "TNF_TREG"

# 每个病人的总T细胞数
patient_total_t <- meta %>%
  group_by(orig.ident) %>%
  summarise(n_Tcell = n(), .groups = "drop")

# 每个病人的目标Treg亚群数
target_treg <- meta %>%
  filter(treg_sub == target_sub) %>%
  group_by(orig.ident) %>%
  summarise(n_target_treg = n(), .groups = "drop")

# 合并
dat_main <- patient_info %>%
  left_join(patient_total_t, by = "orig.ident") %>%
  left_join(target_treg, by = "orig.ident") %>%
  mutate(
    n_target_treg = ifelse(is.na(n_target_treg), 0, n_target_treg),
    prop_target_treg_in_T = n_target_treg / n_Tcell
  )

head(dat_main)
tcell_comp <- meta %>%
  filter(!is.na(T_subtype_refine)) %>%
  group_by(orig.ident, T_subtype_refine) %>%
  summarise(n = n(), .groups = "drop") %>%
  left_join(patient_total_t, by = "orig.ident") %>%
  mutate(prop = n / n_Tcell) %>%
  select(orig.ident, T_subtype_refine, prop) %>%
  pivot_wider(
    names_from = T_subtype_refine,
    values_from = prop,
    values_fill = 0
  )

head(tcell_comp)
treg_only <- meta %>%
  filter(!is.na(treg_sub))

patient_total_treg <- treg_only %>%
  group_by(orig.ident) %>%
  summarise(n_Treg = n(), .groups = "drop")

treg_comp <- treg_only %>%
  group_by(orig.ident, treg_sub) %>%
  summarise(n = n(), .groups = "drop") %>%
  left_join(patient_total_treg, by = "orig.ident") %>%
  mutate(prop = n / n_Treg) %>%
  select(orig.ident, treg_sub, prop) %>%
  pivot_wider(
    names_from = treg_sub,
    values_from = prop,
    values_fill = 0
  )

head(treg_comp)
dat_f <- dat_main %>%
  left_join(tcell_comp, by = "orig.ident")

head(dat_f)
vars_check <- unique(Tcell$T_subtype_refine)
cor_res <- lapply(vars_check, function(v){
  test <- cor.test(dat_f$prop_target_treg_in_T, dat_f[[v]], method = "spearman")
  data.frame(
    celltype = v,
    rho = unname(test$estimate),
    p = test$p.value
  )
}) %>%
  bind_rows() %>%
  mutate(p_adj = p.adjust(p, method = "BH"))

cor_res
library(ggplot2)
library(ggpubr)
library(patchwork)

vars_check <- unique(Tcell$T_subtype_refine)
vars_check <- vars_check[!is.na(vars_check)]

plot_list <- lapply(vars_check, function(v){
  ggplot(dat_f, aes(x = prop_target_treg_in_T, y = .data[[v]])) +
    geom_point(size = 2.5, alpha = 0.85, color = "#2C7FB8") +
    geom_smooth(method = "lm", se = TRUE, color = "#D95F0E", fill = "#FDD0A2") +
    stat_cor(
      method = "spearman",
      label.x.npc = "left",
      label.y.npc = "top",
      size = 4
    ) +
    theme_bw(base_size = 13) +
    labs(
      x = paste0(target_sub, " proportion in T cells"),
      y = paste0(v, " proportion in T cells"),
      title = paste0(target_sub, " vs ", v)
    ) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      axis.title = element_text(face = "bold")
    )
})

wrap_plots(plot_list, ncol = 2)



fit_lm1 <- lm(inflammation_score_num ~ prop_target_treg_in_T, data = dat_f)
summary(fit_lm1)
fit_lm2 <- lm(
  inflammation_score_num ~ prop_target_treg_in_T + CD8_CTL + Trm_CD8 + TH17 + TEM + Activated_TEM_LIKE,
  data = dat_f
)
summary(fit_lm2)
dat_g <- dat_main %>%
  left_join(treg_comp, by = "orig.ident")

head(dat_g)
fit_lm3 <- lm(
  inflammation_score_num ~ TNF_TREG + actived_TREG,
  data = dat_g
)
summary(fit_lm3)
dat_full <- dat_main %>%
  left_join(tcell_comp, by = "orig.ident") %>%
  left_join(treg_comp, by = "orig.ident")

fit_lm4 <- lm(
  inflammation_score_num ~ TNF_TREG + actived_TREG + CD8_CTL + Trm_CD8 + TH17 + TEM + Activated_TEM_LIKE,
  data = dat_full
)
summary(fit_lm4)
library(broom)
library(dplyr)
library(ggplot2)
library(stringr)
library(forcats)

# 整理模型结果
res_lm1 <- tidy(fit_lm1, conf.int = TRUE) %>%
  mutate(model = "Model 1")

res_lm2 <- tidy(fit_lm2, conf.int = TRUE) %>%
  mutate(model = "Model 2")

res_lm3 <- tidy(fit_lm3, conf.int = TRUE) %>%
  mutate(model = "Model 3")

res_lm4 <- tidy(fit_lm4, conf.int = TRUE) %>%
  mutate(model = "Model 4")

res_all <- bind_rows(res_lm1, res_lm2, res_lm3, res_lm4)

# 去掉截距
res_plot <- res_all %>%
  filter(term != "(Intercept)") %>%
  mutate(
    term = recode(term,
                  prop_target_treg_in_T = "TNF_TREG in T cells",
                  TNF_TREG = "TNF_TREG in Tregs",
                  actived_TREG = "Activated TREG",
                  CD8_CTL = "CD8_CTL",
                  Trm_CD8 = "Trm_CD8",
                  TH17 = "TH17",
                  TEM = "TEM",
                  Activated_TEM_LIKE = "Activated_TEM_LIKE"
    )
  )

ggplot(res_plot, aes(x = estimate, y = fct_rev(term), xmin = conf.low, xmax = conf.high)) +
  geom_vline(xintercept = 0, linetype = 2, color = "grey50") +
  geom_errorbarh(height = 0.2, color = "grey40") +
  geom_point(size = 2.8, color = "#2C7FB8") +
  facet_wrap(~model, scales = "free_y") +
  theme_bw(base_size = 13) +
  labs(
    x = "Linear regression coefficient (95% CI)",
    y = NULL,
    title = "Association of TNF_TREG and immune subsets with inflammation score"
  ) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    strip.text = element_text(face = "bold")
  )

tnf_only <- res_all %>%
  filter(term %in% c("prop_target_treg_in_T", "TNF_TREG")) %>%
  mutate(
    exposure = case_when(
      term == "prop_target_treg_in_T" ~ "TNF_TREG in T cells",
      term == "TNF_TREG" ~ "TNF_TREG in Tregs"
    ),
    model = factor(model, levels = c("Model 1", "Model 2", "Model 3", "Model 4"))
  )

ggplot(tnf_only, aes(x = estimate, y = model, xmin = conf.low, xmax = conf.high, color = exposure)) +
  geom_vline(xintercept = 0, linetype = 2, color = "grey50") +
  geom_errorbarh(height = 0.18, position = position_dodge(width = 0.5)) +
  geom_point(size = 3, position = position_dodge(width = 0.5)) +
  theme_bw(base_size = 13) +
  labs(
    x = "Regression coefficient (95% CI)",
    y = NULL,
    color = NULL,
    title = "Effect of TNF_TREG across linear models"
  ) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold")
  )
library(dplyr)
library(ggplot2)

newdata_lm2 <- data.frame(
  prop_target_treg_in_T = seq(
    min(dat_f$prop_target_treg_in_T, na.rm = TRUE),
    max(dat_f$prop_target_treg_in_T, na.rm = TRUE),
    length.out = 100
  ),
  CD8_CTL = median(dat_f$CD8_CTL, na.rm = TRUE),
  Trm_CD8 = median(dat_f$Trm_CD8, na.rm = TRUE),
  TH17 = median(dat_f$TH17, na.rm = TRUE),
  TEM = median(dat_f$TEM, na.rm = TRUE),
  Activated_TEM_LIKE = median(dat_f$Activated_TEM_LIKE, na.rm = TRUE)
)

pred_lm2 <- predict(fit_lm2, newdata = newdata_lm2, interval = "confidence")
pred_df_lm2 <- bind_cols(newdata_lm2, as.data.frame(pred_lm2))

ggplot() +
  geom_point(
    data = dat_f,
    aes(x = prop_target_treg_in_T, y = inflammation_score_num),
    size = 2.5, alpha = 0.7, color = "#2C7FB8"
  ) +
  geom_line(
    data = pred_df_lm2,
    aes(x = prop_target_treg_in_T, y = fit),
    color = "#D95F0E", linewidth = 1.1
  ) +
  geom_ribbon(
    data = pred_df_lm2,
    aes(x = prop_target_treg_in_T, ymin = lwr, ymax = upr),
    alpha = 0.2, fill = "#FDD0A2"
  ) +
  theme_bw(base_size = 13) +
  labs(
    x = "TNF_TREG proportion in T cells",
    y = "Predicted inflammation score",
    title = "Adjusted association between TNF_TREG and inflammation score (Model 2)"
  ) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.title = element_text(face = "bold")
  )


head(IBD_T_CELL_remove_dcs2)
head(TREG_SUBSET)

# 1. 提取 Treg 亚群的关键信息：细胞名 + 细分亚群
treg_anno <- TREG_SUBSET@meta.data[, c("cell_name", "treg_sub")]

# 2. 把亚群信息合并回总 T 细胞群（自动匹配，非 Treg 会是 NA）
# 假设原始大群数据框名为 IBD_T_CELL_remove_dcs2
# Treg亚群数据框名为 TREG_SUBSET，包含列 cell_name 和 treg_sub

library(dplyr)
library(tibble)
treg_anno <- TREG_SUBSET@meta.data[, c("cell_name", "treg_sub")]
library(dplyr)
IBD_T_CELL_remove_dcs2 <- IBD_T_TREG_SUB
IBD_meta <- IBD_T_CELL_remove_dcs2@meta.data
TREG_meta <- TREG_SUBSET@meta.data %>%
  select(cell_name, treg_sub)

# 先去掉旧的 treg_sub，避免生成 .x/.y
IBD_meta <- IBD_meta %>%
  select(-any_of("treg_sub")) %>%
  left_join(TREG_meta, by = "cell_name")

# 非 Treg 补值
IBD_meta$treg_sub[is.na(IBD_meta$treg_sub)] <- "non_Treg"

# 把行名恢复成原来的 Seurat 细胞名
rownames(IBD_meta) <- rownames(IBD_T_CELL_remove_dcs2@meta.data)

IBD_T_CELL_remove_dcs2@meta.data <- IBD_meta

DimPlot(IBD_T_CELL_remove_dcs2, group.by = "treg_sub", label = TRUE)

# 先复制一列最终注释
IBD_T_CELL_remove_dcs2$final_subtype <- IBD_T_CELL_remove_dcs2$T_subtype

# 提取 Treg 子集注释
treg_vec <- TREG_SUBSET@meta.data$treg_sub
names(treg_vec) <- rownames(TREG_SUBSET@meta.data)

# 用 Treg 注释覆盖对应细胞
common_cells <- intersect(rownames(IBD_T_CELL_remove_dcs2@meta.data), names(treg_vec))
IBD_T_CELL_remove_dcs2$final_subtype[common_cells] <- treg_vec[common_cells]

# 画图
DimPlot(IBD_T_CELL_remove_dcs2, group.by = "final_subtype", label = TRUE)

saveRDS(IBD_T_CELL_remove_dcs2,file = 'IBD_T_TREG_SUB.rds')
unique(IBD_T_TREG_SUB$T_subtype)

t_0.8_marker <- FindAllMarkers(IBD_T_CELL_remove_dcs2,
                               group.by = "T_subtype",
                               only.pos = T,min.pct = 0.1,
                               logfc.threshold = 0.5)
write.csv(t_0.8_marker,file = 'IBD_T_CELL_remove_dcs2.CSV')
DotPlot(IBD_T_CELL_remove_dcs2,features = c('CD3E','CD4','CD8A','ITGAE','IL23A','IL7R','TNF'),group.by = 'T_subtype')
TNFRSF4、TNFRSF9（编码4-1BB的基因）和TNFRSF18
DimPlot(IBD_T_CELL_remove_dcs2,group.by = 'T_subtype',label = T)


DotPlot(
  object   = IBD_T_CELL_remove_dcs2,
  group.by = "T_subtype",
  features = marker_genes_list
) + RotatedAxis()

cd4T <- subset(IBD_T_TREG_SUB,T_subtype %>% c("Naive_T","TH17","TFH","Treg") )


###cytotrace
# install.packages("devtools")
# devtools::install_github("digitalcytometry/cytotrace2", subdir = "cytotrace2_r")
remotes::install_github("digitalcytometry/cytotrace2", subdir = "cytotrace2_r")
library(Seurat)
library(CytoTRACE2)
library(ggplot2)
library(dplyr)

# 建议确认默认 assay 是 RNA
DefaultAssay(TREG_SUBSET) <- "RNA"

# 运行 CytoTRACE2
# 人样本用 species = "human"
# 小鼠样本用 species = "mouse"
TREG_SUBSET <- cytotrace2(
  TREG_SUBSET,
  is_seurat = TRUE,
  slot_type = "counts",   # 一般优先用 RNA counts
  species = "human"       # 如果是小鼠改成 "mouse"
)

# 看看新增了哪些 metadata 列
colnames(TREG_SUBSET@meta.data)

# 1) 在原 UMAP 上看连续的 CytoTRACE score
FeaturePlot(
  TREG_SUBSET,
  features = "preKNN_CytoTRACE2_Score",
  reduction = "umap",cols = c('white','red')
)


# 2) 看 potency 分类
DimPlot(
  TREG_SUBSET,
  group.by = "preKNN_CytoTRACE2_Score",
  reduction = "umap",
  label = TRUE
)
DimPlot(
  TREG_SUBSET,
  group.by = "treg_sub",
  reduction = "umap",
  label = TRUE
)
# 3) 继续按你原来的 treg_sub 分组，但比较每个亚群的 CytoTRACE score
VlnPlot(
  TREG_SUBSET,
  features = "CytoTRACE2_Score",
  group.by = "treg_sub",
  pt.size = 0
)

# 4) 也可以做 boxplot
df <- FetchData(TREG_SUBSET, vars = c("treg_sub", "CytoTRACE2_Score"))

ggplot(df, aes(x = treg_sub, y = CytoTRACE2_Score, fill = treg_sub)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.15, size = 0.3, alpha = 0.3) +
  theme_classic() +
  labs(x = "treg_sub", y = "CytoTRACE2 Score")
TREG_SUBSET$CELL_DISEASE <- paste0(TREG_SUBSET$treg_sub,'_',TREG_SUBSET$disease)
DotPlot(TREG_SUBSET,features = c( "MAP2K3", "PMEPA1" ,"RAPGEF6","WNK1",'BATF'),group.by = 'CELL_DISEASE')
library(scKidney)

# 假设您的 Seurat 对象名为 seurat_object

# 指定输出的 H5ad 文件路径
unique(IBD_T_CELL_remove_dcs2$T_subtype)
sce.all.filt <- subset(IBD_T_CELL_remove_dcs2,T_subtype %in% c("Naive_T", "TH17","TFH","Treg" ))
output_h5ad_path <- "CD4T_TREGSUB.h5ad"
setwd("data\\UC")  
# 运行转换函数
getwd()
kidneyH5(sce.all.filt, output_h5ad_path)
DotPlot(TREG_SUBSET,
        features =c('BATF','IL2RA','TNFRSF9','TNFRSF4','TNFRSF18',
                    'FOXP3','IFNG','SELL','CCR7','MAF','FOS','JUN','DUSP1',
                    'IL17A','SLC16A1'),cols = c('white','red'),
        group.by = "treg_sub")+RotatedAxis()
tnf <- subset(TREG_SUBSET,treg_sub=='TNF_TREG')
unique(TREG_SUBSET$treg_sub)
Idents(TREG_SUBSET) <- "disease"
dge <- FindMarkers(
  tnf,
  ident.1 = 'UC',
  ident.2 = 'HC',
  logfc.threshold = 0.5,   # 你原来是 0.5，这里注意参数名是 logfc.threshold
  test.use = "wilcox",
  min.pct = 0.1,
  min.diff.pct = 0.2,
  only.pos = FALSE,        # 关键：要上下调都要，所以设 FALSE
  max.cells.per.ident = 50,
  assay = "RNA"
)
t_marker <- FindAllMarkers(TREG_SUBSET,
                               group.by = "treg_sub",
                               only.pos = T,min.pct = 0.1,
                               logfc.threshold = 0.5)
t_marker %>% filter(cluster=='TNF_TREG'&p_val<0.05&avg_log2FC>1)
tnfr_genes <-  c(
  "TNFRSF18","TNFRSF9","TNFRSF4",
  "LTA","CXCR3","CCL20",
  "IL1R2","IL21R",
  "MALT1","BCL3","PELI1",
  "CFLAR","BCL2L1",
  "MAP2K3","DUSP16",
  "LGALS1"
)

TREG_SUBSET <- subset(TREG_SUBSET,treg_sub=='TNF_TREG')
expr_tnfr <- FetchData(TREG_SUBSET, vars = tnfr_genes)

TREG_SUBSET$TNFR_score <- rowMeans(expr_tnfr, na.rm = TRUE)
genes_use <- c(intersect(res, Treg1$gene_name),'BATF')
genes_use <- genes_use[genes_use %in% rownames(TREG_SUBSET)]

expr <- FetchData(TREG_SUBSET, vars = c("TNFR_score", genes_use))

cor_results <- data.frame(
  gene = genes_use,
  cor = NA,
  pvalue = NA
)

for (g in genes_use) {
  tmp <- cor.test(expr$TNFR_score, expr[[g]], method = "pearson")
  cor_results[cor_results$gene == g, "cor"] <- unname(tmp$estimate)
  cor_results[cor_results$gene == g, "pvalue"] <- tmp$p.value
}

cor_results <- cor_results[order(-cor_results$cor), ]
head(cor_results, 10)
library(ggplot2)

ggplot(cor_results, aes(x = reorder(gene, cor), y = cor)) +
  geom_point(aes(size = -log10(pvalue), color = cor)) +
  coord_flip() +
  theme_bw()

genes_check <- c("MAP2K3","CFLAR","CCL20")

expr <- FetchData(TREG_SUBSET, vars = c("BATF", genes_check))

cor_mat <- cor(expr, method = "spearman")
library(pheatmap)

library(reshape2)
library(ggplot2)

df <- FetchData(TREG_SUBSET, vars = c("BATF", "CFLAR"))
df <- df[df$BATF  & df$CFLAR, ]
cor_mat <- cor(df, method = "spearman", use = "pairwise.complete.obs")
cor_mat
ggplot(df, aes(BATF, CFLAR)) +
  geom_point(alpha = 0.4) +
  geom_smooth(method = "lm", color = "red") +
  theme_classic() +
  ggtitle("BATF vs MAP2K3")


df <- FetchData(
  TREG_SUBSET,
  vars = c("CFLAR", "BATF", "orig.ident")
)

# 每个病人聚合成一个值
patient_df <- df %>%
  group_by(orig.ident) %>%
  summarise(
    CFLAR = mean(CFLAR, na.rm = TRUE),
    BATF  = mean(BATF,  na.rm = TRUE),
    ncell = n()
  ) %>%
  ungroup()
test <- cor.test(patient_df$CFLAR, patient_df$BATF, method = "pearson")

test
ggplot(patient_df, aes(x = CFLAR, y = BATF)) +
  geom_point( alpha = 0.8) +
  geom_smooth(method = "lm", se = FALSE) +
  annotate(
    "text",
    x = Inf, y = Inf,
    label = paste0(
      "Pearson r = ", round(test$estimate, 3),
      "\np = ", signif(test$p.value, 3)
    ),
    hjust = 1.1, vjust = 1.5,
    size = 5
  ) +
  theme_classic() +
  labs(
    x = "Mean FOXP3 per patient",
    y = "Mean BATF per patient",
    title = "Patient-level correlation",
    size = "Cell number"
  )
library(reshape2)
library(ggplot2)
ggplot(df, aes(BATF, CFLAR)) +
  geom_point(alpha = 0.3) +
  geom_smooth(method = "lm", color = "red") +
  stat_cor(method = "spearman") +   # 需要 ggpubr
  theme_classic()



df <- FetchData(
  TREG_SUBSET,
  vars = c("FOXP3", "BATF", "orig.ident")
)

# 每个病人聚合成一个值
patient_df <- df %>%
  group_by(orig.ident) %>%
  summarise(
    FOXP3 = mean(FOXP3, na.rm = TRUE),
    BATF  = mean(BATF,  na.rm = TRUE),
    ncell = n()
  ) %>%
  ungroup()
library(reshape2)
library(ggplot2)
ggplot(patient_df, aes(BATF, FOXP3)) +
  geom_point(alpha = 0.3) +
  geom_smooth(method = "lm", color = "red") +
  stat_cor(method = "spearman") +   # 需要 ggpubr
  theme_classic()
DimPlot(
  ,
  group.by = "final_subtype",
  reduction = "umap",
  label = TRUE
)
unique(IBD_T_TREG_SUB$final_subtype)
unique(TREG_SUBSET $treg_sub)
