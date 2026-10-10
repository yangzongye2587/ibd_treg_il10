# ========= 仅需运行一次：加载包 + 全局函数 =========
setwd("data\\bulkil10\\过程文件")


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
    
    # symbol转ENTREZ
    entrez_all <- mapIds(
      org.Mm.eg.db,
      keys = dt_clean$gene.symbol,
      column = "ENTREZID",
      keytype = "SYMBOL",
      multiVals = "first"
    )
    
    # 筛选有效ID，同步保留对应symbol（解决names为空行数不匹配）
    valid_index <- !is.na(entrez_all)
    valid_entrez <- unique(entrez_all[valid_index])
    valid_symbol <- dt_clean$gene.symbol[valid_index]
    
    if (length(valid_entrez) < 3) {
      warning(paste0(prefix, " 可转换ENTREZ基因不足3个，跳过富集"))
      return()
    }
    cat("有效ENTREZ基因数：", length(valid_entrez), "\n")
    
    # 构建正确长度匹配的映射表
    id_map <- data.frame(
      geneID = as.character(valid_entrez),
      gene.symbol = valid_symbol,
      stringsAsFactors = FALSE
    )
    
    # KEGG富集（标准阈值）
    kegg_res <- enrichKEGG(
      gene = valid_entrez,
      organism = "mmu",
      pvalueCutoff = 0.95,
      qvalueCutoff = 0.95
    )
    
    # GO三大分类富集
    go_bp <- enrichGO(gene = valid_entrez, OrgDb = org.Mm.eg.db, ont = "BP", pAdjustMethod = "BH", pvalueCutoff = 0.05, qvalueCutoff = 0.2)
    go_mf <- enrichGO(gene = valid_entrez, OrgDb = org.Mm.eg.db, ont = "MF", pAdjustMethod = "BH", pvalueCutoff = 0.05, qvalueCutoff = 0.2)
    go_cc <- enrichGO(gene = valid_entrez, OrgDb = org.Mm.eg.db, ont = "CC", pAdjustMethod = "BH", pvalueCutoff = 0.05, qvalueCutoff = 0.2)
    
    # 创建输出文件夹
    out_dir <- paste0(prefix, "_FJ")
    if (!dir.exists(out_dir)) dir.create(out_dir)
    
    # 导出KEGG结果+绘图
    if (nrow(kegg_res) > 0) {
      kegg_df <- as.data.frame(kegg_res)
      kegg_df$gene.symbol <- sapply(kegg_df$geneID, function(x) tryCatch(get_gene_symbols(x, id_map), error = function(e) NA))
      write.csv(kegg_df, file.path(out_dir, "KEGG.csv"), row.names = F)
      
      pdf(file.path(out_dir, "KEGG_dotplot.pdf"), width = 8, height = 6)
      print(dotplot(kegg_res, title = paste0(prefix, " KEGG富集")))
      dev.off()
    } else {
      warning(paste0(prefix, " KEGG无显著富集条目"))
    }
    
    # GO BP
    if (nrow(go_bp) > 0) {
      bp_df <- as.data.frame(go_bp)
      bp_df$gene.symbol <- sapply(bp_df$geneID, function(x) tryCatch(get_gene_symbols(x, id_map), error = function(e) NA))
      write.csv(bp_df, file.path(out_dir, "GO_BP.csv"), row.names = F)
      
      pdf(file.path(out_dir, "GO_BP_dotplot.pdf"), width = 8, height = 6)
      print(dotplot(go_bp, title = paste0(prefix, " GO-BP富集")))
      dev.off()
    } else {
      warning(paste0(prefix, " GO-BP无显著富集条目"))
    }
    
    # GO MF
    if (nrow(go_mf) > 0) {
      mf_df <- as.data.frame(go_mf)
      mf_df$gene.symbol <- sapply(mf_df$geneID, function(x) tryCatch(get_gene_symbols(x, id_map), error = function(e) NA))
      write.csv(mf_df, file.path(out_dir, "GO_MF.csv"), row.names = F)
      
      pdf(file.path(out_dir, "GO_MF_dotplot.pdf"), width = 8, height = 6)
      print(dotplot(go_mf, title = paste0(prefix, " GO-MF富集")))
      dev.off()
    } else {
      warning(paste0(prefix, " GO-MF无显著富集条目"))
    }
    
    # GO CC
    if (nrow(go_cc) > 0) {
      cc_df <- as.data.frame(go_cc)
      cc_df$gene.symbol <- sapply(cc_df$geneID, function(x) tryCatch(get_gene_symbols(x, id_map), error = function(e) NA))
      write.csv(cc_df, file.path(out_dir, "GO_CC.csv"), row.names = F)
      
      pdf(file.path(out_dir, "GO_CC_dotplot.pdf"), width = 8, height = 6)
      print(dotplot(go_cc, title = paste0(prefix, " GO-CC富集")))
      dev.off()
    } else {
      warning(paste0(prefix, " GO-CC无显著富集条目"))
    }
    
    cat(prefix, " 处理完成，输出目录：", out_dir, "\n")
  }) # local结束，当前循环所有临时变量彻底销毁
  
  # 强制垃圾回收，释放内存，杜绝下一轮循环被上一轮残留干扰
  gc(verbose = FALSE)
  cat("<<<<<<", prefix, "内存清理完毕，进入下一组\n\n")
}

cat("======== 全部文件富集流程执行完毕 ========")
