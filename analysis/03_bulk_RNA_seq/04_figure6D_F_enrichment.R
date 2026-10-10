# ========= 仅需运行一次：加载包 + 全局函数 =========
setwd("E:/Y")
library(clusterProfiler)
library(org.Mm.eg.db)
library(pathview)
library(dplyr)
library(ggplot2)
# 安全基因ID匹配函数（全局，无依赖污染）
get_gene_symbols <- function(gene_id_str, map_df) {
  if (is.na(gene_id_str) || gene_id_str == "" || gene_id_str == "/") {
    return(NA)
  }
  id_vec <- unlist(strsplit(as.character(gene_id_str), "/"))
  id_vec <- id_vec[id_vec != ""]
  if (length(id_vec) == 0) return(NA)
  
  hit_sym <- map_df$gene.symbol[map_df$geneID %in% id_vec]
  if (length(hit_sym) == 0) NA else paste(unique(hit_sym), collapse = "/")
}
# 待处理文件前缀列表
file_prefix <- c(
  "TN_vs_IBD_UP", "TN_vs_IBD_DOWN", "TN_vs_IBD_ALL",
  "IN_vs_IBD_UP", "IN_vs_IBD_DOWN", "IN_vs_IBD_ALL",
  "TN_vs_IN_UP", "TN_vs_IN_DOWN", "TN_vs_IN_ALL"
)
# ========= 单层for循环，单次处理完全隔离 =========
for (prefix in file_prefix) {
  cat("\n=============================================\n")
  cat("开始处理：", prefix, "\n")
  
  # 单次循环所有变量封闭在local内部，循环结束自动销毁，无交叉干扰
  local({
    csv_path <- paste0(prefix, ".csv")
    if (!file.exists(csv_path)) {
      warning(paste("文件缺失，跳过：", csv_path))
      return()
    }
    
    # 读取并清洗Excel空白脏行
    dt_raw <- read.csv(csv_path, stringsAsFactors = FALSE)
    dt_clean <- dt_raw[!(is.na(dt_raw$gene.symbol) | dt_raw$gene.symbol == ""), ]
    if (nrow(dt_clean) < 3) {
      warning(paste0(prefix, " 有效基因不足3个，跳过富集"))
      return()
    }
    
    # -------------------------- 基因ID转换逻辑 --------------------------
    # symbol转ENTREZID（与原逻辑完全一致）
    gene_sym <- dt_clean$gene.symbol
    gene_map <- bitr(gene_sym,
                     fromType = "SYMBOL",
                     toType   = "ENTREZID",
                     OrgDb    = org.Mm.eg.db)
    genes_entrez <- unique(gene_map$ENTREZID)
    
    if (length(genes_entrez) < 3) {
      warning(paste0(prefix, " 可转换ENTREZ基因不足3个，跳过富集"))
      return()
    }
    cat("有效ENTREZ基因数：", length(genes_entrez), "\n")
    
    # 构建ID映射表（用于后续geneID转symbol）
    id_map <- data.frame(
      geneID = as.character(gene_map$ENTREZID),
      gene.symbol = gene_map$SYMBOL,
      stringsAsFactors = FALSE
    )
    
    # -------------------------- GO富集 + 新增TOP8筛选+ID截断 --------------------------
    # GO富集（ONT=ALL，与原逻辑一致）
    eG <- enrichGO(
      gene         = genes_entrez,
      OrgDb        = org.Mm.eg.db,
      pvalueCutoff = 0.05,
      ont          = "ALL",
      readable     = TRUE
    )
    
    # ========= 新增：按Count降序取TOP8 + 基因ID截断（最多前10个） =========
    eG_df <- as.data.frame(eG) %>%
      filter(!is.na(p.adjust)) %>%        # 先去NA
      arrange(-Count) %>%               # 按Count降序排序
      slice_head(n = 8) %>%              # 取TOP8条目
      mutate(
        # 基因ID截断：按/分割后取前10个，超出的截断
        geneID_short = sapply(
          strsplit(geneID, "/"),
          function(x) paste(head(x, 10), collapse = "/")
        )
      )
    
    # -------------------------- KEGG富集 + 新增TOP8筛选+ID截断 --------------------------
    # KEGG富集（参数与原逻辑完全一致）
    KEGG_ids <- mapIds(x       = org.Mm.eg.db,
                       keys    = gene_sym,
                       keytype = "SYMBOL",
                       column  = "ENTREZID")
    KEGG_ids <- na.omit(KEGG_ids)
    
    kk <- enrichKEGG(
      gene = KEGG_ids,
      organism = "mmu",
      keyType = "kegg",
      pAdjustMethod = "fdr",
      pvalueCutoff = 0.05
    )
    
    # ========= 新增：按Count降序取TOP8 + 基因ID截断（最多前10个） =========
    ekegg_readable <- setReadable(kk,
                                  OrgDb  = org.Mm.eg.db,
                                  keyType = "ENTREZID")
    kegg_filtered <- as.data.frame(ekegg_readable) %>%
      filter(pvalue < 0.05, Count > 5) %>%
      filter(category != 'Human Diseases') %>%
      arrange(-Count) %>%  # 按Count降序排序
      slice_head(n = 8) %>% # 取TOP8条目
      mutate(
        # 基因ID截断：按/分割后取前10个，超出的截断
        geneID_short = sapply(
          strsplit(geneID, "/"),
          function(x) paste(head(x, 10), collapse = "/")
        )
      )
    
    # -------------------------- 输出文件夹创建 --------------------------
    out_dir <- paste0(prefix, "_FJ")
    if (!dir.exists(out_dir)) dir.create(out_dir)
    
    # -------------------------- 导出GO结果 + 可视化（适配截断后的ID） --------------------------
    if (nrow(eG_df) > 0) {
      # 导出GO结果（同时保留原始geneID和截断后的geneID_short）
      write.csv(eG_df, file.path(out_dir, "GO_ALL.csv"), row.names = F)
      
      # GO可视化（完全适配原定制主题，仅替换为截断后的ID）
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
      # x轴从0起
      p1 <- p + scale_x_continuous(expand = c(0, 0))
      # 画Description文本
      p2 <- p1 +
        theme(axis.text.y = element_blank()) +
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
      # 叠加截断后的geneID信息（核心修改：用geneID_short替代原始geneID）
      p3 <- p2 +
        geom_text(
          data = eG_df,
          aes(
            x = 0.1,
            y = Description,
            label = geneID_short,
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
          title = paste0(prefix, "_GO富集_TOP8"),
          x     = "Count",
          y     = "Description"
        )
      # 应用主题和配色
      p5 <- p4 + mytheme +
        scale_fill_manual(values = mycol) +
        scale_color_manual(values = mycol)
      
      # 导出GO可视化PDF
      pdf(file.path(out_dir, "GO_ALL_barplot_TOP8.pdf"), width = 12, height = 10)
      print(p5)
      dev.off()
    } else {
      warning(paste0(prefix, " GO-ALL无符合筛选条件的富集条目"))
    }
    
    # -------------------------- 导出KEGG结果 + 可视化（适配截断后的ID） --------------------------
    if (nrow(kegg_filtered) > 0) {
      # 导出KEGG结果（同时保留原始geneID和截断后的geneID_short）
      write.csv(kegg_filtered, file.path(out_dir, "KEGG.csv"), row.names = F)
      
      # KEGG可视化（完全适配原定制主题，仅替换为截断后的ID）
      mytheme <- theme(
        plot.title   = element_text(size = 20, face = "bold"),
        axis.title   = element_text(size = 18),
        axis.text    = element_text(size = 16),
        axis.ticks.y = element_blank()
      )
      mycol <- c("#6fbedd", "#d9543d", "#54B946", "#F89E1B", "#BC5428FF")
      
      p <- ggplot() +
        geom_bar(
          data = kegg_filtered,
          aes(x = Count, y = reorder(Description, Count), fill = category),
          width = 0.5,
          stat = "identity"
        ) +
        theme_classic()
      # x轴从0起
      p1 <- p + scale_x_continuous(expand = c(0, 0))
      # 画Description文本
      p2 <- p1 +
        theme(axis.text.y = element_blank()) +
        geom_text(
          data = kegg_filtered,
          aes(
            x = 0.1,
            y = Description,
            label = Description
          ),
          size = 8,
          hjust = 0
        )
      # 叠加截断后的geneID信息（核心修改：用geneID_short替代原始geneID）
      p3 <- p2 +
        geom_text(
          data = kegg_filtered,
          aes(
            x = 0.1,
            y = Description,
            label = geneID_short,
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
          title = paste0(prefix, "_KEGG富集_TOP8"),
          x     = "Count",
          y     = "Description"
        )
      # 应用主题和配色
      p5 <- p4 + mytheme +
        scale_fill_manual(values = mycol) +
        scale_color_manual(values = mycol)
      
      # 导出KEGG可视化PDF
      pdf(file.path(out_dir, "KEGG_barplot_TOP8.pdf"), width = 12, height = 10)
      print(p5)
      dev.off()
    } else {
      warning(paste0(prefix, " KEGG无符合筛选条件的富集条目"))
    }
    
    cat(prefix, " 处理完成，输出目录：", out_dir, "\n")
  }) # local结束，当前循环所有临时变量彻底销毁
  
  # 强制垃圾回收，释放内存，杜绝下一轮循环被上一轮残留干扰
  gc(verbose = FALSE)
  cat("<<<<<<", prefix, "内存清理完毕，进入下一组\n\n")
}
cat("======== 全部文件富集流程执行完毕 ========")
