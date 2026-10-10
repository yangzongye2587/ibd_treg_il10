cellRatioPlot(object = TREG_SUBSET,
              sample.name = "disease",
              celltype.name = "treg_sub")
unique(TREG_SUBSET$treg_sub)
####marker的kegg/go####
TREG_SUBSET <- readRDS("data/TREG_SUBSET.RDS")
aaa <- FindMarkers(object = TREG_SUBSET,ident.1 = 'TNF_TREG',ident.2 = 'actived_TREG',
                   group.by = 'treg_sub')
aaa <- aaa %>% filter(p_val<0.05&avg_log2FC>1)


suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(stringr)
  library(org.Hs.eg.db)
  library(clusterProfiler)
  library(org.Mm.eg.db)   # 小鼠
  # library(org.Hs.eg.db) # 人类的话用这个
})
aaa$gene <- rownames(aaa)
gene_sym <- aaa$gene       # SYMBOL 列

gene_map <- bitr(gene_sym,
                 fromType = "SYMBOL",
                 toType   = "ENTREZID",
                 OrgDb    = org.Hs.eg.db)

genes_entrez <- unique(gene_map$ENTREZID)

eG <- enrichGO(
  gene         = genes_entrez,
  OrgDb        = org.Hs.eg.db,
  pvalueCutoff = 0.05,
  ont          = "ALL",   # 或 "BP" / "CC" / "MF"
  readable     = TRUE
)
go <- as.data.frame(eG) 
eG_df <- as.data.frame(eG) %>%
  filter(!is.na(p.adjust)) %>%        # **** 先去 NA
  arrange(-Count) %>%               # **** 按 p.adjust 从小到大
  slice_head(n = 8)                  # **** 只保留前 15 条

mytheme <- theme(
  plot.title   = element_text(size = 20, face = "bold"),
  axis.title   = element_text(size = 18),
  axis.text    = element_text(size = 16),
  axis.ticks.y = element_blank()
)

mycol <- c("#6fbedd", "#d9543d", "#54B946", "#F89E1B", "#BC5428FF")
p <- ggplot() +
  geom_bar(
    data = eG_df,
    aes(x = Count, y = reorder(Description, Count), fill = ONTOLOGY),
    width = 0.5,
    stat = "identity"
  ) +
  theme_classic()

# x 轴从 0 起
p1 <- p + scale_x_continuous(expand = c(0, 0))

# 先画 Description 文本
p2 <- p1 +
  theme(axis.text.y = element_blank()) +  # 去掉 y 轴标签
  geom_text(
    data = eG_df,
    aes(
      x = 0.1,
      y = Description,
      label = Description
    ),
    size = 8,
    hjust = 0
  )

# 叠加 geneID 信息
p3 <- p2 +
  geom_text(
    data = eG_df %>%
      mutate(
        geneID = sapply(geneID, function(gene_str) {
          # 按 '/' 分割基因ID字符串
          gene_list <- strsplit(gene_str, "/")[[1]]
          # 选取前10个基因，合并为一个字符串
          paste(head(gene_list, 10), collapse = "/")
        })
      ),
    aes(
      x = 0.1,
      y = Description,
      label = geneID,
      color = ONTOLOGY
    ),
    size = 6,
    fontface = "italic",
    hjust = 0,
    vjust = 2.7
  )

# 标题和坐标轴标签
p4 <- p3 +
  labs(
    title = 'GO_TNFR_TREG_MARKER_GENE',         # group 名作为主标题
    x     = "Count",
    y     = "Description"
  )

# 应用主题和配色
p5 <- p4 + mytheme +
  scale_fill_manual(values = mycol) +
  scale_color_manual(values = mycol)


print(p5)



KEGG_ids <- mapIds(x       = org.Hs.eg.db,
                   keys    = gene_sym,
                   keytype = "SYMBOL",
                   column  = "ENTREZID")
KEGG_ids <- na.omit(KEGG_ids)



kk <- enrichKEGG(
  gene = KEGG_ids,          # 或者用上面 ncbi-geneid 的写法
  organism = "hsa",
  keyType = "kegg"
)



ekegg_readable <- setReadable(kk,
                              OrgDb  = org.Hs.eg.db,
                              keyType = "ENTREZID")


down <- as.data.frame(ekegg_readable)
down <- down %>% filter(pvalue<0.05,Count>3)
down <- down %>% filter(category !='Human Diseases')
mytheme <- theme(
  plot.title   = element_text(size = 20, face = "bold"),
  axis.title   = element_text(size = 18),
  axis.text    = element_text(size = 16),
  axis.ticks.y = element_blank()
)

mycol <- c("#6fbedd", "#d9543d", "#54B946", "#F89E1B", "#BC5428FF")
p <- ggplot() +
  geom_bar(
    data = down,
    aes(x = Count, y = reorder(Description, Count), fill = category),
    width = 0.5,
    stat = "identity"
  ) +
  theme_classic()

# x 轴从 0 起
p1 <- p + scale_x_continuous(expand = c(0, 0))

# 先画 Description 文本
p2 <- p1 +
  theme(axis.text.y = element_blank()) +  # 去掉 y 轴标签
  geom_text(
    data = down,
    aes(
      x = 0.1,
      y = Description,
      label = Description
    ),
    size = 8,
    hjust = 0
  )

# 叠加 geneID 信息
p3 <- p2 +
  geom_text(
    data = down %>%
      mutate(
        geneID = sapply(geneID, function(gene_str) {
          # 按 '/' 分割基因ID字符串
          gene_list <- strsplit(gene_str, "/")[[1]]
          # 选取前10个基因，合并为一个字符串
          paste(head(gene_list, 10), collapse = "/")
        })
      ),
    aes(
      x = 0.1,
      y = Description,
      label = geneID,
      color = category
    ),
    size = 6,
    fontface = "italic",
    hjust = 0,
    vjust = 2.7
  )

# 标题和坐标轴标签
p4 <- p3 +
  labs(
    title = 'KEGG_TNFR_TREG_VS_ACTIVED_TREG',         # group 名作为主标题
    x     = "Count",
    y     = "Description"
  )

# 应用主题和配色
p5 <- p4 + mytheme +
  scale_fill_manual(values = mycol) +
  scale_color_manual(values = mycol)


print(p5)

####df的kegg/go####
TREG_SUBSET <- readRDS("data/TREG_SUBSET.RDS")
obj <- subset(TREG_SUBSET,treg_sub=='TNF_TREG')
unique(TREG_SUBSET$disease)
aaa <- FindMarkers(object = obj,ident.1 = 'UC',ident.2 ='HC' ,group.by = 'disease')
aaa <- aaa %>% filter(p_val<0.05&avg_log2FC>1)


suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(stringr)
  library(org.Hs.eg.db)
  library(clusterProfiler)
  library(org.Mm.eg.db)   # 小鼠
  # library(org.Hs.eg.db) # 人类的话用这个
})
aaa$gene <- rownames(aaa)
gene_sym <- aaa$gene       # SYMBOL 列

gene_map <- bitr(gene_sym,
                 fromType = "SYMBOL",
                 toType   = "ENTREZID",
                 OrgDb    = org.Hs.eg.db)

genes_entrez <- unique(gene_map$ENTREZID)

eG <- enrichGO(
  gene         = genes_entrez,
  OrgDb        = org.Hs.eg.db,
  pvalueCutoff = 0.05,
  ont          = "ALL",   # 或 "BP" / "CC" / "MF"
  readable     = TRUE
)

eG_df <- as.data.frame(eG) %>%
  filter(!is.na(p.adjust)) %>%        # **** 先去 NA
  arrange(-Count) %>%               # **** 按 p.adjust 从小到大
  slice_head(n = 10)                  # **** 只保留前 15 条

mytheme <- theme(
  plot.title   = element_text(size = 20, face = "bold"),
  axis.title   = element_text(size = 18),
  axis.text    = element_text(size = 16),
  axis.ticks.y = element_blank()
)

mycol <- c("#6fbedd", "#d9543d", "#54B946", "#F89E1B", "#BC5428FF")
p <- ggplot() +
  geom_bar(
    data = eG_df,
    aes(x = Count, y = reorder(Description, Count), fill = ONTOLOGY),
    width = 0.5,
    stat = "identity"
  ) +
  theme_classic()

# x 轴从 0 起
p1 <- p + scale_x_continuous(expand = c(0, 0))

# 先画 Description 文本
p2 <- p1 +
  theme(axis.text.y = element_blank()) +  # 去掉 y 轴标签
  geom_text(
    data = eG_df,
    aes(
      x = 0.1,
      y = Description,
      label = Description
    ),
    size = 8,
    hjust = 0
  )

# 叠加 geneID 信息
p3 <- p2 +
  geom_text(
    data = eG_df %>%
      mutate(
        geneID = sapply(geneID, function(gene_str) {
          # 按 '/' 分割基因ID字符串
          gene_list <- strsplit(gene_str, "/")[[1]]
          # 选取前10个基因，合并为一个字符串
          paste(head(gene_list, 10), collapse = "/")
        })
      ),
    aes(
      x = 0.1,
      y = Description,
      label = geneID,
      color = ONTOLOGY
    ),
    size = 6,
    fontface = "italic",
    hjust = 0,
    vjust = 2.7
  )

# 标题和坐标轴标签
p4 <- p3 +
  labs(
    title = 'GO_TNFR_TREG_df_GENE',         # group 名作为主标题
    x     = "Count",
    y     = "Description"
  )

# 应用主题和配色
p5 <- p4 + mytheme +
  scale_fill_manual(values = mycol) +
  scale_color_manual(values = mycol)


print(p5)



KEGG_ids <- mapIds(x       = org.Hs.eg.db,
                   keys    = gene_sym,
                   keytype = "SYMBOL",
                   column  = "ENTREZID")
KEGG_ids <- na.omit(KEGG_ids)



kk <- enrichKEGG(
  gene = KEGG_ids,          # 或者用上面 ncbi-geneid 的写法
  organism = "hsa",
  keyType = "kegg"
)



ekegg_readable <- setReadable(kk,
                              OrgDb  = org.Hs.eg.db,
                              keyType = "ENTREZID")


down <- as.data.frame(ekegg_readable)
down <- down %>% filter(pvalue<0.05,Count>3)
down <- down %>% filter(category !='Human Diseases')
mytheme <- theme(
  plot.title   = element_text(size = 20, face = "bold"),
  axis.title   = element_text(size = 18),
  axis.text    = element_text(size = 16),
  axis.ticks.y = element_blank()
)

mycol <- c("#6fbedd", "#d9543d", "#54B946", "#F89E1B", "#BC5428FF")
p <- ggplot() +
  geom_bar(
    data = down,
    aes(x = Count, y = reorder(Description, Count), fill = category),
    width = 0.5,
    stat = "identity"
  ) +
  theme_classic()

# x 轴从 0 起
p1 <- p + scale_x_continuous(expand = c(0, 0))

# 先画 Description 文本
p2 <- p1 +
  theme(axis.text.y = element_blank()) +  # 去掉 y 轴标签
  geom_text(
    data = down,
    aes(
      x = 0.1,
      y = Description,
      label = Description
    ),
    size = 8,
    hjust = 0
  )

# 叠加 geneID 信息
p3 <- p2 +
  geom_text(
    data = down %>%
      mutate(
        geneID = sapply(geneID, function(gene_str) {
          # 按 '/' 分割基因ID字符串
          gene_list <- strsplit(gene_str, "/")[[1]]
          # 选取前10个基因，合并为一个字符串
          paste(head(gene_list, 10), collapse = "/")
        })
      ),
    aes(
      x = 0.1,
      y = Description,
      label = geneID,
      color = category
    ),
    size = 6,
    fontface = "italic",
    hjust = 0,
    vjust = 2.7
  )

# 标题和坐标轴标签
p4 <- p3 +
  labs(
    title = 'KEGG_TNFR_TREG_df_GENE',         # group 名作为主标题
    x     = "Count",
    y     = "Description"
  )

# 应用主题和配色
p5 <- p4 + mytheme +
  scale_fill_manual(values = mycol) +
  scale_color_manual(values = mycol)


print(p5)

####df_2CELL的kegg/go####
TREG_SUBSET <- readRDS("data/TREG_SUBSET.RDS")
obj <- subset(TREG_SUBSET,treg_sub=='TNF_TREG')
unique(TREG_SUBSET$treg_sub)
aaa <- FindMarkers(obj = TREG_SUBSET,ident.1 = 'TNF_TREG',ident.2 ='actived_TREG' ,group.by = 'treg_sub')
aaa <- aaa %>% filter(p_val<0.05&avg_log2FC>1)


suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(stringr)
  library(org.Hs.eg.db)
  library(clusterProfiler)
  library(org.Mm.eg.db)   # 小鼠
  # library(org.Hs.eg.db) # 人类的话用这个
})
aaa$gene <- rownames(aaa)
gene_sym <- aaa$gene       # SYMBOL 列

gene_map <- bitr(gene_sym,
                 fromType = "SYMBOL",
                 toType   = "ENTREZID",
                 OrgDb    = org.Hs.eg.db)

genes_entrez <- unique(gene_map$ENTREZID)

eG <- enrichGO(
  gene         = genes_entrez,
  OrgDb        = org.Hs.eg.db,
  pvalueCutoff = 0.05,
  ont          = "ALL",   # 或 "BP" / "CC" / "MF"
  readable     = TRUE
)

eG_df <- as.data.frame(eG) %>%
  filter(!is.na(p.adjust)) %>%        # **** 先去 NA
  arrange(-Count) %>%               # **** 按 p.adjust 从小到大
  slice_head(n = 8)                  # **** 只保留前 15 条

mytheme <- theme(
  plot.title   = element_text(size = 20, face = "bold"),
  axis.title   = element_text(size = 18),
  axis.text    = element_text(size = 16),
  axis.ticks.y = element_blank()
)

mycol <- c("#6fbedd", "#d9543d", "#54B946", "#F89E1B", "#BC5428FF")
p <- ggplot() +
  geom_bar(
    data = eG_df,
    aes(x = Count, y = reorder(Description, Count), fill = ONTOLOGY),
    width = 0.5,
    stat = "identity"
  ) +
  theme_classic()

# x 轴从 0 起
p1 <- p + scale_x_continuous(expand = c(0, 0))

# 先画 Description 文本
p2 <- p1 +
  theme(axis.text.y = element_blank()) +  # 去掉 y 轴标签
  geom_text(
    data = eG_df,
    aes(
      x = 0.1,
      y = Description,
      label = Description
    ),
    size = 8,
    hjust = 0
  )

# 叠加 geneID 信息
p3 <- p2 +
  geom_text(
    data = eG_df %>%
      mutate(
        geneID = sapply(geneID, function(gene_str) {
          # 按 '/' 分割基因ID字符串
          gene_list <- strsplit(gene_str, "/")[[1]]
          # 选取前10个基因，合并为一个字符串
          paste(head(gene_list, 10), collapse = "/")
        })
      ),
    aes(
      x = 0.1,
      y = Description,
      label = geneID,
      color = ONTOLOGY
    ),
    size = 6,
    fontface = "italic",
    hjust = 0,
    vjust = 2.7
  )

# 标题和坐标轴标签
p4 <- p3 +
  labs(
    title = 'GO_TNFR_TREG_VS_ACTIVED_TREG',         # group 名作为主标题
    x     = "Count",
    y     = "Description"
  )

# 应用主题和配色
p5 <- p4 + mytheme +
  scale_fill_manual(values = mycol) +
  scale_color_manual(values = mycol)


print(p5)

