library(readr)
library(tidyverse)
setwd("data")  
mono_drivers_TREG <- read_csv("mono_drivers_TREG.csv")
treg_top_hubs <- read_csv("UC\\treg_top_hubs.csv")
drivers <- mono_drivers_TREG %>% filter(Treg_pval<0.05&Treg_corr>0.2)
treg_wgcna <- treg_top_hubs %>% filter(module=='Treg1'&kME_Treg1>0.2)
inter <- intersect(treg_wgcna$gene_name,drivers$gene_name)

DotPlot(IBD_T_CELL_remove_dcs2,
        features = c(treg_wgcna$gene_name),group.by = 'CELL_DISE')+RotatedAxis()












####GSE75214_RAW.tar####
library(affy)
library(simpleaffy)
cel_dir = "data/bulk_validation\\GSE75214_RAW" 
raw <- ReadAffy(celfile.path = cel_dir)  # 读入全部CEL
# 后续常规：RMA 标准化
eset <- rma(raw)
expr <- exprs(eset)  # 表达矩阵

library(AnnotationDbi)
library(hugene10sttranscriptcluster.db)

probe_ids <- rownames(expr)   # 你这里是 7892501 这种

# 1) 先自动找出哪个 keytype 能匹配你的行名（避免再报 key 无效）
kt <- keytypes(hugene10sttranscriptcluster.db)
test <- probe_ids[1:min(200, length(probe_ids))]
score <- sapply(kt, function(k){
  sum(test %in% keys(hugene10sttranscriptcluster.db, keytype = k))
})
best_keytype <- names(sort(score, decreasing = TRUE))[1]
best_keytype
sort(score, decreasing=TRUE)[1:10]
# 2) 用最匹配的 keytype 做映射
map <- select(
  hugene10sttranscriptcluster.db,
  keys     = probe_ids,
  keytype  = best_keytype,
  columns  = c("SYMBOL", "ENTREZID")
)
# 3) 对齐回表达矩阵，并按基因聚合（同一基因多个探针取平均）
map <- map[!is.na(map$SYMBOL) & map$SYMBOL != "", ]
expr2 <- expr[match(map[[best_keytype]], rownames(expr)), , drop=FALSE]
gene  <- map$SYMBOL
expr_gene_sum <- rowsum(expr2, group = gene)
expr_gene <- expr_gene_sum / as.vector(table(gene))

expr_gene_GSE75214 <- expr_gene
####转基因名
expr_gene_GSE75214 <- expr_gene  # 你已经有了

old_names <- colnames(expr_gene_GSE75214)

# 从列名里提取 GSM号
gsm <- sub(".*(GSM[0-9]+).*", "\\1", old_names)

# 检查：有没有提取失败（失败的话会变成原字符串）
bad <- !grepl("^GSM[0-9]+$", gsm)
if (any(bad)) {
  print(old_names[bad][1:min(10, sum(bad))])
  stop("有些列名里没找到 GSM 号，请检查列名格式。")
}

# 检查：有没有重复 GSM（一般不该重复）
if (any(duplicated(gsm))) {
  dup <- gsm[duplicated(gsm)]
  print(unique(dup))
  stop("提取到重复 GSM，说明列名里可能包含相同 GSM 或文件重复。")
}

# 替换列名
colnames(expr_gene_GSE75214) <- gsm
head(colnames(expr_gene_GSE75214))

head(expr_gene_GSE75214)
gse <- getGEO("GSE75214", GSEMatrix = TRUE)
eset_geo <- gse[[1]]

clinical <- pData(eset_geo)
clinical$sample <- rownames(clinical)
###验证diff
setdiff(clinical$sample,colnames(expr_gene_GSE75214))
###提取想要的clinical

clinical1 <- clinical %>% filter(`disease:ch1` %in% c('ulcerative colitis','control'))
table(clinical1$`disease:ch1`)
gsm_selct <- clinical1$sample
expr_gene_GSE75214_filter <- expr_gene_GSE75214[,gsm_selct]
clinical1
clinical1$phenotype <- NA
clinical1$phenotype[clinical1$`disease:ch1` == "ulcerative colitis"] <- "uc"
clinical1$phenotype[clinical1$`disease:ch1` == "control"] <- "hc"



####GSE87473_raw####
library(affy)
cel_dir = "data/bulk_validation\\GSE87473_raw" 
raw <- ReadAffy(celfile.path = cel_dir)  # 读入全部CEL
# 后续常规：RMA 标准化
eset <- rma(raw)
expr <- exprs(eset)  # 表达矩阵

probe_ids <- rownames(expr)   # 你这里是 7892501 这种


library(AnnotationDbi)

expr <- exprs(eset)
probe_pm <- rownames(expr)

# 1) 去掉 _PM_，把 1007_PM_s_at -> 1007_s_at
probe <- sub("_PM_", "_", probe_pm)

# 2) 备选注释包（按常见优先级：U133 Plus2 / U133A / PrimeView）
cand <- list(
  hgu133plus2 = hgu133plus2.db
)

# 3) 选一个“命中率最高”的注释包
score_one <- function(db){
  kt <- keytypes(db)
  if (!("PROBEID" %in% kt)) return(0)
  test <- probe[1:min(200, length(probe))]
  sum(test %in% keys(db, keytype="PROBEID"))
}

scores <- sapply(cand, score_one)
scores
best_name <- names(which.max(scores))
best_db <- cand[[best_name]]
best_name

# 4) 映射 probe set id -> gene symbol
map <- AnnotationDbi::select(
  best_db,
  keys     = probe,
  keytype  = "PROBEID",
  columns  = c("SYMBOL", "ENTREZID")
)

# 5) 对齐回 expr，并按 SYMBOL 聚合（同一基因多个探针取平均）
map <- map[!is.na(map$SYMBOL) & map$SYMBOL != "", ]
expr2 <- expr[match(sub("_PM_", "_", map$PROBEID), probe), , drop=FALSE]
gene  <- map$SYMBOL

expr_gene_sum <- rowsum(expr2, group = gene)
expr_gene <- expr_gene_sum / as.vector(table(gene))

dim(expr_gene)
head(rownames(expr_gene))

GSE87473_raw_data_expr_gene<- expr_gene 
####转基因名
GSE87473_raw_data_expr_gene <- expr_gene  # 你已经有了

old_names <- colnames(GSE87473_raw_data_expr_gene)

# 从列名里提取 GSM号
gsm <- sub(".*(GSM[0-9]+).*", "\\1", old_names)

# 检查：有没有提取失败（失败的话会变成原字符串）
bad <- !grepl("^GSM[0-9]+$", gsm)
if (any(bad)) {
  print(old_names[bad][1:min(10, sum(bad))])
  stop("有些列名里没找到 GSM 号，请检查列名格式。")
}

# 检查：有没有重复 GSM（一般不该重复）
if (any(duplicated(gsm))) {
  dup <- gsm[duplicated(gsm)]
  print(unique(dup))
  stop("提取到重复 GSM，说明列名里可能包含相同 GSM 或文件重复。")
}

# 替换列名
colnames(GSE87473_raw_data_expr_gene) <- gsm
head(colnames(GSE87473_raw_data_expr_gene))

head(GSE87473_raw_data_expr_gene)
gse <- getGEO("GSE87473", GSEMatrix = TRUE, AnnotGPL = FALSE)
eset_geo <- gse[[1]]

clinical <- pData(eset_geo)
clinical$sample <- rownames(clinical)
###验证diff
setdiff(clinical$sample,colnames(GSE87473_raw_data_expr_gene))
###提取想要的clinical
table(clinical$`disease:ch1`)
clinical_GSE87473 <- clinical

# 新建 phenotype 列
clinical_GSE87473$phenotype <- NA
clinical_GSE87473$phenotype[clinical_GSE87473$`disease:ch1` == "Ulcerative Colitis (UC)"] <- "uc"
clinical_GSE87473$phenotype[clinical_GSE87473$`disease:ch1` == "Normal"] <- "hc"
gsm_selct <- clinical_GSE87473$sample
GSE87473_raw_data_expr_gene_filter <- GSE87473_raw_data_expr_gene[,gsm_selct]


####train_data_test_data####
####train_data
test.x= expr_gene_GSE75214_filter
test.y = clinical1[,c("sample","phenotype")]

###test_data
train.x = GSE87473_raw_data_expr_gene_filter 
colnames(clinical_GSE87473)
train.y= clinical_GSE87473[,c("sample","phenotype")]
####gene_hub
treg_wgcna <- treg_top_hubs %>% filter(module=='Treg1'&kME_Treg1>0.25)
gene_hub <- treg_wgcna$gene_name
gene_hub
head(train.x)
head(train.y)
head(test.x)
head(test.y)



########训练####

##只有lasso
genes <- gene_hub
g_use <- intersect(genes, intersect(rownames(train.x), rownames(test.x)))

X_train <- t(train.x[g_use, , drop=FALSE])
X_test  <- t(test.x[g_use, , drop=FALSE])

y_train <- factor(train.y$phenotype, levels=c("hc","uc"))
y_test  <- factor(test.y$phenotype,  levels=c("hc","uc"))

# 去 NA
keep_tr <- !is.na(y_train); keep_te <- !is.na(y_test)
X_train <- X_train[keep_tr, , drop=FALSE]; y_train <- y_train[keep_tr]
X_test  <- X_test[keep_te, , drop=FALSE]; y_test  <- y_test[keep_te]
library(glmnet)

set.seed(1)
ytr <- ifelse(y_train=="uc", 1, 0)

# 标准化（只用训练集参数）
mu <- apply(X_train,2,mean)
sdv <- apply(X_train,2,sd); sdv[sdv==0] <- 1
Xtr <- scale(X_train, center=mu, scale=sdv)
Xte <- scale(X_test,  center=mu, scale=sdv)

cv <- cv.glmnet(as.matrix(Xtr), ytr, family="binomial", alpha=1, nfolds=10, type.measure="auc")
fit_lasso <- glmnet(as.matrix(Xtr), ytr, family="binomial", alpha=1, lambda=cv$lambda.min)

sel <- rownames(as.matrix(coef(fit_lasso)))
sel <- sel[as.numeric(coef(fit_lasso)) != 0]
sel <- setdiff(sel, "(Intercept)")
sel

####直接用这些算法（并做重要性/排序）

library(ranger)
library(pROC)

Xtr_use <- Xtr  # 或 Xtr[, sel]
Xte_use <- Xte  # 或 Xte[, sel]

df_tr <- data.frame(y=y_train, Xtr_use)
df_te <- data.frame(y=y_test,  Xte_use)

set.seed(1)
rf <- ranger(y ~ ., data=df_tr, probability=TRUE, num.trees=1000,
             importance="permutation", classification=TRUE)

pred <- predict(rf, data=df_te)$predictions[, "uc"]
auc <- as.numeric(roc(ifelse(y_test=="uc",1,0), pred, quiet=TRUE)$auc)
auc

# 基因重要性（越大越重要）
imp <- sort(rf$variable.importance, decreasing=TRUE)
imp <- as.data.frame.Date(imp)





#####SVM（RBF核）+ permutation importance（推荐）




library(e1071)
library(pROC)

set.seed(1)
svmfit <- svm(x=as.data.frame(Xtr), y=y_train, kernel="radial", probability=TRUE)

pred <- attr(predict(svmfit, as.data.frame(Xte), probability=TRUE), "probabilities")[, "uc"]
auc <- as.numeric(roc(ifelse(y_test=="uc",1,0), pred, quiet=TRUE)$auc)
auc









library(xgboost)
library(pROC)

dtrain <- xgb.DMatrix(data=as.matrix(Xtr), label=ifelse(y_train=="uc",1,0))
dtest  <- xgb.DMatrix(data=as.matrix(Xte), label=ifelse(y_test=="uc",1,0))

set.seed(1)
bst <- xgb.train(
  data=dtrain, nrounds=200,
  params=list(
    objective="binary:logistic",
    eval_metric="auc",
    max_depth=3,
    eta=0.05,
    subsample=0.8,
    colsample_bytree=0.8
  ),
  verbose=0
)

pred <- predict(bst, dtest)
auc <- as.numeric(roc(ifelse(y_test=="uc",1,0), pred, quiet=TRUE)$auc)
auc

# 重要性排序
imp <- xgb.importance(model=bst)
imp
library(ggplot2)
library(data.table)
library(pROC)
# imp: xgb.importance 的结果（data.table / data.frame）
imp_dt <- as.data.table(imp)
imp_dt <- imp_dt[order(-Gain)]
topN <- 20  # 想看前多少个自己改

imp_top <- imp_dt[1:min(topN, .N)]
imp_top[, Feature := factor(Feature, levels=rev(Feature))]

p1 <- ggplot(imp_top, aes(x=Feature, y=Gain)) +
  geom_col() +
  coord_flip() +
  labs(title=paste0("XGBoost Feature Importance (Gain) - Top ", topN),
       x=NULL, y="Gain")

print(p1)
imp_dt[, cumGain := cumsum(Gain)]
imp_dt[, cumPct := cumGain / sum(Gain)]
imp_dt[, rank := .I]

# 找到 90/95/99% 的截断点
n90 <- imp_dt[cumPct >= 0.90, rank][1]
n95 <- imp_dt[cumPct >= 0.95, rank][1]
n99 <- imp_dt[cumPct >= 0.99, rank][1]

p2 <- ggplot(imp_dt, aes(x=rank, y=cumPct)) +
  geom_line() +
  geom_point(size=1) +
  geom_hline(yintercept=c(0.90,0.95,0.99), linetype=2) +
  geom_vline(xintercept=c(n90,n95,n99), linetype=2) +
  scale_y_continuous(limits=c(0,1)) +
  labs(title="Cumulative Gain (XGBoost)",
       x="Top-k features (ranked by Gain)", y="Cumulative proportion of total Gain")

print(p2)

cat("Top-k for 90%:", n90, "\n")
cat("Top-k for 95%:", n95, "\n")
cat("Top-k for 99%:", n99, "\n")
# y_test: factor(hc, uc)
# pred_te: numeric 概率
roc_te <- roc(ifelse(y_test=="uc", 1, 0), pred_te, quiet=TRUE)
auc_te <- as.numeric(roc_te$auc)

roc_df <- data.frame(
  fpr = 1 - roc_te$specificities,
  tpr = roc_te$sensitivities
)

p3 <- ggplot(roc_df, aes(x=fpr, y=tpr)) +
  geom_line() +
  geom_abline(linetype=2) +
  coord_equal() +
  labs(title=paste0("External Test ROC (AUC=", round(auc_te, 3), ")"),
       x="False Positive Rate (1 - Specificity)",
       y="True Positive Rate (Sensitivity)")

print(p3)



####训练画####
library(glmnet)
library(ranger)
library(e1071)
library(xgboost)
library(pROC)
library(ggplot2)
library(data.table)
## ===== LASSO 图 =====

# 1) CV 曲线（glmnet自带）
plot(cv)
abline(v = log(c(cv$lambda.min, cv$lambda.1se)), lty = 2, col = c("red", "blue"))
title("LASSO CV Curve")

# 2) 系数路径图（glmnet自带）
plot(fit_lasso, xvar = "lambda", label = TRUE)
title("LASSO Coefficient Path")

# 3) 非零系数条形图
coef_lasso <- coef(fit_lasso)
coef_df <- data.frame(
  Feature = rownames(coef_lasso),
  Coef = as.numeric(coef_lasso)
)
coef_df <- subset(coef_df, Feature != "(Intercept)" & Coef != 0)

coef_df$Feature <- factor(coef_df$Feature, levels = coef_df$Feature[order(coef_df$Coef)])

p_lasso_coef <- ggplot(coef_df, aes(x = Feature, y = Coef)) +
  geom_col() +
  coord_flip() +
  labs(title = "LASSO Non-zero Coefficients", x = NULL, y = "Coefficient")

print(p_lasso_coef)

## ===== Random Forest =====

Xtr_use <- Xtr   # 或 Xtr[, sel, drop=FALSE]
Xte_use <- Xte   # 或 Xte[, sel, drop=FALSE]

df_tr <- data.frame(y = y_train, Xtr_use)
df_te <- data.frame(y = y_test,  Xte_use)

set.seed(1)
rf <- ranger(
  y ~ ., data = df_tr,
  probability = TRUE,
  num.trees = 1000,
  importance = "permutation",
  classification = TRUE
)

pred_rf <- predict(rf, data = df_te)$predictions[, "uc"]
auc_rf <- as.numeric(roc(ifelse(y_test == "uc", 1, 0), pred_rf, quiet = TRUE)$auc)
auc_rf
imp_rf <- sort(rf$variable.importance, decreasing = TRUE)
imp_rf_df <- data.frame(
  Feature = names(imp_rf),
  Importance = as.numeric(imp_rf)
)

topN <- 20
imp_rf_top <- head(imp_rf_df, topN)
imp_rf_top$Feature <- factor(imp_rf_top$Feature, levels = rev(imp_rf_top$Feature))

p_rf_imp <- ggplot(imp_rf_top, aes(x = Feature, y = Importance)) +
  geom_col() +
  coord_flip() +
  labs(title = paste0("Random Forest Variable Importance (Top ", topN, ")"),
       x = NULL, y = "Permutation Importance")

print(p_rf_imp)
roc_rf <- roc(ifelse(y_test == "uc", 1, 0), pred_rf, quiet = TRUE)

roc_rf_df <- data.frame(
  fpr = 1 - roc_rf$specificities,
  tpr = roc_rf$sensitivities
)

# 加载必要的包
library(ggplot2)
library(ggthemes) # 提供更多美观的主题
library(scales)   # 用于刻度格式化

# 美化后的 ROC 曲线代码
p_rf_roc <- ggplot(roc_rf_df, aes(x = fpr, y = tpr)) +
  # 核心线条 - 自定义颜色、粗细和透明度
  geom_line(color = "#2E86AB", linewidth = 1.2, alpha = 0.8) +
  # 参考线（对角线）- 调整样式和颜色
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", 
              color = "#A23B72", linewidth = 0.8, alpha = 0.6) +
  # 等比例坐标
  coord_equal(ratio = 1, xlim = c(0, 1), ylim = c(0, 1)) +
  # 标题和坐标轴标签（优化字体大小和样式）
  labs(
    title = paste0("Random Forest ROC Curve"),
    subtitle = paste0("AUC = ", round(auc_rf, 3)),
    x = "False Positive Rate (FPR)",
    y = "True Positive Rate (TPR)",
    caption = "Data: Random Forest Classification Results"
  ) +
  # 主题设置（核心美化部分）
  theme_minimal() +
  theme(
    # 标题样式
    plot.title = element_text(size = 16, face = "bold", color = "#2C3E50", hjust = 0.5),
    plot.subtitle = element_text(size = 14, color = "#34495E", hjust = 0.5),
    plot.caption = element_text(size = 10, color = "#7F8C8D", hjust = 1),
    # 坐标轴样式
    axis.title.x = element_text(size = 12, color = "#2C3E50", margin = margin(t = 10)),
    axis.title.y = element_text(size = 12, color = "#2C3E50", margin = margin(r = 10)),
    axis.text = element_text(size = 10, color = "#34495E"),
    axis.line = element_line(color = "#BDC3C7", linewidth = 0.5),
    # 网格线
    panel.grid.major = element_line(color = "#ECF0F1", linewidth = 0.3),
    panel.grid.minor = element_blank(),
    # 背景
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA),
    # 边距
    plot.margin = margin(20, 20, 20, 20)
  ) +
  # 坐标轴刻度优化
  scale_x_continuous(breaks = seq(0, 1, 0.1), limits = c(0, 1), expand = c(0, 0)) +
  scale_y_continuous(breaks = seq(0, 1, 0.1), limits = c(0, 1), expand = c(0, 0))

# 如果需要添加 AUC 数值标注（可选）
p_rf_roc <- p_rf_roc +
  annotate("text", x = 0.7, y = 0.2, 
           label = paste0("AUC = ", round(auc_rf, 3)),
           size = 5, color = "#2E86AB", fontface = "bold")

# 显示图表
print(p_rf_roc)

# 保存高清图表（可选）
ggsave("rf_roc_curve.png", p_rf_roc, width = 8, height = 7, dpi = 300, bg = "white")

print(p_rf_roc)
## ===== SVM (RBF) =====



pred_svm <- attr(
  predict(svmfit, as.data.frame(Xte), probability = TRUE),
  "probabilities"
)[, "uc"]

auc_svm <- as.numeric(roc(ifelse(y_test == "uc", 1, 0), pred_svm, quiet = TRUE)$auc)
auc_svm
roc_svm <- roc(ifelse(y_test == "uc", 1, 0), pred_svm, quiet = TRUE)

roc_svm_df <- data.frame(
  fpr = 1 - roc_svm$specificities,
  tpr = roc_svm$sensitivities
)

p_svm_roc <- ggplot(roc_svm_df, aes(x = fpr, y = tpr)) +
  geom_line() +
  geom_abline(linetype = 2) +
  coord_equal() +
  labs(title = paste0("SVM ROC (AUC = ", round(auc_svm, 3), ")"),
       x = "False Positive Rate",
       y = "True Positive Rate")

print(p_svm_roc)
## permutation importance for SVM
set.seed(1)

base_auc <- auc_svm
features <- colnames(Xte)
perm_imp <- numeric(length(features))
names(perm_imp) <- features

for (j in seq_along(features)) {
  Xte_perm <- Xte
  Xte_perm[, j] <- sample(Xte_perm[, j])
  
  pred_perm <- attr(
    predict(svmfit, as.data.frame(Xte_perm), probability = TRUE),
    "probabilities"
  )[, "uc"]
  
  auc_perm <- as.numeric(roc(ifelse(y_test == "uc", 1, 0), pred_perm, quiet = TRUE)$auc)
  perm_imp[j] <- base_auc - auc_perm
}

imp_svm_df <- data.frame(
  Feature = names(perm_imp),
  Importance = perm_imp
)
imp_svm_df <- imp_svm_df[order(-imp_svm_df$Importance), ]

topN <- 20
imp_svm_top <- head(imp_svm_df, topN)
imp_svm_top$Feature <- factor(imp_svm_top$Feature, levels = rev(imp_svm_top$Feature))

p_svm_imp <- ggplot(imp_svm_top, aes(x = Feature, y = Importance)) +
  geom_col() +
  coord_flip() +
  labs(title = paste0("SVM Permutation Importance (Top ", topN, ")"),
       x = NULL, y = "AUC Drop After Permutation")

print(p_svm_imp)
## ===== XGBoost =====

dtest  <- xgb.DMatrix(data = as.matrix(Xte), label = ifelse(y_test == "uc", 1, 0))

pred_xgb <- predict(bst, dtest)
auc_xgb <- as.numeric(roc(ifelse(y_test == "uc", 1, 0), pred_xgb, quiet = TRUE)$auc)
auc_xgb
imp_xgb <- xgb.importance(model = bst)
imp_dt <- as.data.table(imp_xgb)
imp_dt <- imp_dt[order(-Gain)]

topN <- 20
imp_top <- imp_dt[1:min(topN, .N)]
imp_top[, Feature := factor(Feature, levels = rev(Feature))]

p_xgb_imp <- ggplot(imp_top, aes(x = Feature, y = Gain)) +
  geom_col() +
  coord_flip() +
  labs(title = paste0("XGBoost Feature Importance (Gain) - Top ", topN),
       x = NULL, y = "Gain")

print(p_xgb_imp)

imp_dt[, cumGain := cumsum(Gain)]
imp_dt[, cumPct := cumGain / sum(Gain)]
imp_dt[, rank := .I]

n90 <- imp_dt[cumPct >= 0.90, rank][1]
n95 <- imp_dt[cumPct >= 0.95, rank][1]
n99 <- imp_dt[cumPct >= 0.99, rank][1]

p_xgb_cum <- ggplot(imp_dt, aes(x = rank, y = cumPct)) +
  geom_line() +
  geom_point(size = 1) +
  geom_hline(yintercept = c(0.90, 0.95, 0.99), linetype = 2) +
  geom_vline(xintercept = c(n90, n95, n99), linetype = 2) +
  scale_y_continuous(limits = c(0, 1)) +
  labs(title = "Cumulative Gain (XGBoost)",
       x = "Top-k features (ranked by Gain)",
       y = "Cumulative proportion of total Gain")

print(p_xgb_cum)

cat("Top-k for 90%:", n90, "\n")
cat("Top-k for 95%:", n95, "\n")
cat("Top-k for 99%:", n99, "\n")
roc_xgb <- roc(ifelse(y_test == "uc", 1, 0), pred_xgb, quiet = TRUE)
auc_xgb <- as.numeric(roc_xgb$auc)

roc_xgb_df <- data.frame(
  fpr = 1 - roc_xgb$specificities,
  tpr = roc_xgb$sensitivities
)

p_xgb_roc <- ggplot(roc_xgb_df, aes(x = fpr, y = tpr)) +
  geom_line() +
  geom_abline(linetype = 2) +
  coord_equal() +
  labs(title = paste0("External Test ROC (AUC = ", round(auc_xgb, 3), ")"),
       x = "False Positive Rate (1 - Specificity)",
       y = "True Positive Rate (Sensitivity)")

print(p_xgb_roc)
## ===== Combined ROC =====

roc_lasso_pred <- as.numeric(predict(fit_lasso, newx = as.matrix(Xte), type = "response"))
roc_lasso_obj <- roc(ifelse(y_test == "uc", 1, 0), roc_lasso_pred, quiet = TRUE)
auc_lasso <- as.numeric(roc_lasso_obj$auc)

roc_rf_obj <- roc(ifelse(y_test == "uc", 1, 0), pred_rf, quiet = TRUE)
roc_svm_obj <- roc(ifelse(y_test == "uc", 1, 0), pred_svm, quiet = TRUE)
roc_xgb_obj <- roc(ifelse(y_test == "uc", 1, 0), pred_xgb, quiet = TRUE)

roc_all_df <- rbind(
  data.frame(
    fpr = 1 - roc_lasso_obj$specificities,
    tpr = roc_lasso_obj$sensitivities,
    Model = paste0("LASSO (AUC=", round(as.numeric(roc_lasso_obj$auc), 3), ")")
  ),
  data.frame(
    fpr = 1 - roc_rf_obj$specificities,
    tpr = roc_rf_obj$sensitivities,
    Model = paste0("RF (AUC=", round(as.numeric(roc_rf_obj$auc), 3), ")")
  ),
  data.frame(
    fpr = 1 - roc_svm_obj$specificities,
    tpr = roc_svm_obj$sensitivities,
    Model = paste0("SVM (AUC=", round(as.numeric(roc_svm_obj$auc), 3), ")")
  ),
  data.frame(
    fpr = 1 - roc_xgb_obj$specificities,
    tpr = roc_xgb_obj$sensitivities,
    Model = paste0("XGBoost (AUC=", round(as.numeric(roc_xgb_obj$auc), 3), ")")
  )
)

p_all_roc <- ggplot(roc_all_df, aes(x = fpr, y = tpr, color = Model)) +
  geom_line(linewidth = 1) +
  geom_abline(linetype = 2) +
  coord_equal() +
  labs(title = "ROC Comparison on External Test Set",
       x = "False Positive Rate",
       y = "True Positive Rate")

print(p_all_roc)



####画图精修
library(ggplot2)
library(grid)

theme_cns <- function(base_size = 14, base_family = "sans") {
  theme_classic(base_size = base_size, base_family = base_family) +
    theme(
      plot.title = element_text(size = base_size + 1, face = "bold", hjust = 0.5),
      axis.title = element_text(size = base_size, face = "bold", color = "black"),
      axis.text = element_text(size = base_size - 1, color = "black"),
      axis.line = element_line(color = "black", linewidth = 0.8),
      axis.ticks = element_line(color = "black", linewidth = 0.6),
      axis.ticks.length = unit(0.18, "cm"),
      legend.title = element_text(size = base_size - 1, face = "bold"),
      legend.text = element_text(size = base_size - 2),
      legend.position = "right",
      panel.border = element_blank(),
      plot.margin = margin(10, 12, 10, 10)
    )
}

coef_lasso <- coef(fit_lasso)
coef_df <- data.frame(
  Feature = rownames(coef_lasso),
  Coef = as.numeric(coef_lasso)
)
coef_df <- subset(coef_df, Feature != "(Intercept)" & Coef != 0)
coef_df <- coef_df[order(coef_df$Coef), ]
coef_df$Feature <- factor(coef_df$Feature, levels = coef_df$Feature)

p_lasso_coef <- ggplot(coef_df, aes(x = Feature, y = Coef)) +
  geom_col(width = 0.72, fill = "#4C78A8", color = "black", linewidth = 0.25) +
  coord_flip() +
  labs(
    title = "LASSO Selected Features",
    x = NULL,
    y = "Coefficient"
  ) +
  theme_cns()

print(p_lasso_coef)
imp_rf <- sort(rf$variable.importance, decreasing = TRUE)
imp_rf_df <- data.frame(
  Feature = names(imp_rf),
  Importance = as.numeric(imp_rf)
)

topN <- 20
imp_rf_top <- head(imp_rf_df, topN)
imp_rf_top <- imp_rf_top[order(imp_rf_top$Importance), ]
imp_rf_top$Feature <- factor(imp_rf_top$Feature, levels = imp_rf_top$Feature)

p_rf_imp <- ggplot(imp_rf_top, aes(x = Feature, y = Importance)) +
  geom_col(width = 0.72, fill = "#4C78A8", color = "black", linewidth = 0.25) +
  coord_flip() +
  labs(
    title = paste0("Random Forest Importance (Top ", topN, ")"),
    x = NULL,
    y = "Permutation Importance"
  ) +
  theme_cns()

print(p_rf_imp)
imp_svm_df <- imp_svm_df[order(imp_svm_df$Importance, decreasing = FALSE), ]
topN <- 20
imp_svm_top <- tail(imp_svm_df, topN)
imp_svm_top <- imp_svm_top[order(imp_svm_top$Importance), ]
imp_svm_top$Feature <- factor(imp_svm_top$Feature, levels = imp_svm_top$Feature)

p_svm_imp <- ggplot(imp_svm_top, aes(x = Feature, y = Importance)) +
  geom_col(width = 0.72, fill = "#4C78A8", color = "black", linewidth = 0.25) +
  coord_flip() +
  labs(
    title = paste0("SVM Permutation Importance (Top ", topN, ")"),
    x = NULL,
    y = "AUC Decrease"
  ) +
  theme_cns()

print(p_svm_imp)
imp_xgb <- xgb.importance(model = bst)
imp_dt <- as.data.table(imp_xgb)
imp_dt <- imp_dt[order(-Gain)]

topN <- 20
imp_top <- imp_dt[1:min(topN, .N)]
imp_top <- imp_top[order(Gain)]
imp_top[, Feature := factor(Feature, levels = Feature)]

p_xgb_imp <- ggplot(imp_top, aes(x = Feature, y = Gain)) +
  geom_col(width = 0.72, fill = "#4C78A8", color = "black", linewidth = 0.25) +
  coord_flip() +
  labs(
    title = paste0("XGBoost Feature Importance (Top ", topN, ")"),
    x = NULL,
    y = "Gain"
  ) +
  theme_cns()

print(p_xgb_imp)



plot_roc_cns <- function(roc_obj, auc_value, title_text, line_color = "#D62728") {
  roc_df <- data.frame(
    fpr = 1 - roc_obj$specificities,
    tpr = roc_obj$sensitivities
  )
  
  ggplot(roc_df, aes(x = fpr, y = tpr)) +
    geom_line(linewidth = 1.3, color = line_color) +
    geom_abline(intercept = 0, slope = 1, linetype = 2, color = "grey55", linewidth = 0.7) +
    annotate(
      "text", x = 0.62, y = 0.12,
      label = paste0("AUC = ", sprintf("%.3f", auc_value)),
      size = 5, fontface = "bold"
    ) +
    coord_equal() +
    scale_x_continuous(limits = c(0, 1), expand = c(0, 0)) +
    scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) +
    labs(
      title = title_text,
      x = "False Positive Rate",
      y = "True Positive Rate"
    ) +
    theme_cns()
}
pred_lasso <- as.numeric(predict(fit_lasso, newx = as.matrix(Xte), type = "response"))
roc_lasso <- roc(ifelse(y_test == "uc", 1, 0), pred_lasso, quiet = TRUE)
auc_lasso <- as.numeric(roc_lasso$auc)

p_lasso_roc <- plot_roc_cns(
  roc_obj = roc_lasso,
  auc_value = auc_lasso,
  title_text = "LASSO ROC",
  line_color = "#D62728"
)

print(p_lasso_roc)


roc_rf <- roc(ifelse(y_test == "uc", 1, 0), pred_rf, quiet = TRUE)
auc_rf <- as.numeric(roc_rf$auc)

p_rf_roc <- plot_roc_cns(
  roc_obj = roc_rf,
  auc_value = auc_rf,
  title_text = "Random Forest ROC",
  line_color = "#D62728"
)

print(p_rf_roc)
roc_svm <- roc(ifelse(y_test == "uc", 1, 0), pred_svm, quiet = TRUE)
auc_svm <- as.numeric(roc_svm$auc)

p_svm_roc <- plot_roc_cns(
  roc_obj = roc_svm,
  auc_value = auc_svm,
  title_text = "SVM ROC",
  line_color = "#D62728"
)

print(p_svm_roc)
roc_xgb <- roc(ifelse(y_test == "uc", 1, 0), pred_xgb, quiet = TRUE)
auc_xgb <- as.numeric(roc_xgb$auc)

p_xgb_roc <- plot_roc_cns(
  roc_obj = roc_xgb,
  auc_value = auc_xgb,
  title_text = "XGBoost ROC",
  line_color = "#D62728"
)

print(p_xgb_roc)
roc_rf_obj    <- roc(ifelse(y_test == "uc", 1, 0), pred_rf, quiet = TRUE)
roc_svm_obj   <- roc(ifelse(y_test == "uc", 1, 0), pred_svm, quiet = TRUE)
roc_xgb_obj   <- roc(ifelse(y_test == "uc", 1, 0), pred_xgb, quiet = TRUE)
roc_lasso_obj <- roc(ifelse(y_test == "uc", 1, 0), pred_lasso, quiet = TRUE)

roc_all_df <- rbind(
  data.frame(
    fpr = 1 - roc_lasso_obj$specificities,
    tpr = roc_lasso_obj$sensitivities,
    Model = paste0("LASSO (AUC=", sprintf("%.3f", as.numeric(roc_lasso_obj$auc)), ")")
  ),
  data.frame(
    fpr = 1 - roc_rf_obj$specificities,
    tpr = roc_rf_obj$sensitivities,
    Model = paste0("RF (AUC=", sprintf("%.3f", as.numeric(roc_rf_obj$auc)), ")")
  ),
  data.frame(
    fpr = 1 - roc_svm_obj$specificities,
    tpr = roc_svm_obj$sensitivities,
    Model = paste0("SVM (AUC=", sprintf("%.3f", as.numeric(roc_svm_obj$auc)), ")")
  ),
  data.frame(
    fpr = 1 - roc_xgb_obj$specificities,
    tpr = roc_xgb_obj$sensitivities,
    Model = paste0("XGBoost (AUC=", sprintf("%.3f", as.numeric(roc_xgb_obj$auc)), ")")
  )
)

p_all_roc <- ggplot(roc_all_df, aes(x = fpr, y = tpr, color = Model)) +
  geom_line(linewidth = 1.15) +
  geom_abline(intercept = 0, slope = 1, linetype = 2, color = "grey55", linewidth = 0.7) +
  coord_equal() +
  scale_x_continuous(limits = c(0, 1), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) +
  scale_color_manual(values = c(
    "LASSO (AUC=0.000)" = "#1F77B4",
    "RF (AUC=0.000)" = "#D62728",
    "SVM (AUC=0.000)" = "#2CA02C",
    "XGBoost (AUC=0.000)" = "#9467BD"
  )) +
  labs(
    title = "ROC Curves of Different Models",
    x = "False Positive Rate",
    y = "True Positive Rate",
    color = NULL
  ) +
  theme_cns() +
  theme(
    legend.position = c(0.68, 0.25),
    legend.background = element_blank()
  )

print(p_all_roc)
model_levels <- unique(roc_all_df$Model)
model_cols <- c("#1F77B4", "#D62728", "#2CA02C", "#9467BD")
names(model_cols) <- model_levels

p_all_roc <- ggplot(roc_all_df, aes(x = fpr, y = tpr, color = Model)) +
  geom_line(linewidth = 1.15) +
  geom_abline(intercept = 0, slope = 1, linetype = 2, color = "grey55", linewidth = 0.7) +
  coord_equal() +
  scale_x_continuous(limits = c(0, 1), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) +
  scale_color_manual(values = model_cols) +
  labs(
    title = "ROC Curves of Different Models",
    x = "False Positive Rate",
    y = "True Positive Rate",
    color = NULL
  ) +
  theme_cns() +
  theme(
    legend.position = c(0.68, 0.25),
    legend.background = element_blank()
  )

print(p_all_roc)


remotes::install_github("souravc83/fastAdaboost")

remotes::install_github("l-magnificence/Mime")
