
library(DESeq2)
setwd("E:/Y/")

exprSet <- read.table("TN_vs_IBD.txt", head = T)
head(exprSet)

group <- factor(c("TN","TN","TN",
                  "control","control","control"),
                levels = c("control","TN"))

Data <- data.frame(row.names = colnames(exprSet), group = group)
head(Data)

dds <- DESeqDataSetFromMatrix(countData = exprSet,
                              colData = Data,
                              design = ~ group)

dds <- DESeq(dds)

sizeFactors(dds)

# ==== 修改： results() ====
res <- results(dds)

## 归一化数据保存
resdata <- merge(as.data.frame(res), 
                 as.data.frame(counts(dds, normalized = TRUE)),
                 by = "row.names", sort = FALSE)
head(resdata)
write.csv(resdata, file = "DESeq2_TN_vs_IBD.csv", row.names = FALSE)

#################################################



library(DESeq2)
setwd("E:/Y/")

exprSet <- read.table("IN_vs_IBD.txt", head = T)
head(exprSet)

group <- factor(c("IN","IN","IN",
                  "control","control","control"),
                levels = c("control","IN"))

Data <- data.frame(row.names = colnames(exprSet), group = group)
head(Data)

dds <- DESeqDataSetFromMatrix(countData = exprSet,
                              colData = Data,
                              design = ~ group)

dds <- DESeq(dds)

sizeFactors(dds)

# ==== 修改： results() ====
res <- results(dds)

## 归一化数据保存
resdata <- merge(as.data.frame(res), 
                 as.data.frame(counts(dds, normalized = TRUE)),
                 by = "row.names", sort = FALSE)
head(resdata)
write.csv(resdata, file = "DESeq2_IN_vs_IBD.csv", row.names = FALSE)



##########################################################


library(DESeq2)
setwd("E:/Y/")

exprSet <- read.table("TN_vs_IN.txt", head = T)
head(exprSet)

group <- factor(c("TN","TN","TN",
                  "control","control","control"),
                levels = c("control","TN"))

Data <- data.frame(row.names = colnames(exprSet), group = group)
head(Data)

dds <- DESeqDataSetFromMatrix(countData = exprSet,
                              colData = Data,
                              design = ~ group)

dds <- DESeq(dds)

sizeFactors(dds)

# ==== 修改： results() ====
res <- results(dds)

## 归一化数据保存
resdata <- merge(as.data.frame(res), 
                 as.data.frame(counts(dds, normalized = TRUE)),
                 by = "row.names", sort = FALSE)
head(resdata)
write.csv(resdata, file = "DESeq2_TN_vs_IN.csv", row.names = FALSE)

