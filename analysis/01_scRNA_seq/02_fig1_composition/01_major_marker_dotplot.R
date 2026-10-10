library(Seurat)
library(dplyr)
library(ggplot2)
library(grid)
obj <- major_annotation_ibd_subtype_tregback
# =========================
# 1. 每个细胞类型扩充到 3 个 marker
# =========================
library(Seurat)
library(dplyr)
library(ggplot2)
library(grid)

# =========================
# 1. 细胞类型顺序
# =========================
celltype_order <- c(
  "B",
  "endothelial_cells",
  "Epithelial",
  "glia",
  "Mast",
  "MYLOID",
  "neuron",
  "other",
  "pericyte",
  "plasma",
  "Stromal",
  "T/NK"
)

obj$cell_subtype <- factor(
  obj$cell_subtype,
  levels = celltype_order
)

# =========================
# 2. 每个细胞类型对应 3 个 marker
# =========================
marker_genes_list <- list(
  "B" = c("CD19", "MS4A1", "CD79A"),
  
  "endothelial_cells" = c("PECAM1", "VWF", "CDH5"),
  
  "Epithelial" = c("EPCAM", "KRT8", "KRT18"),
  
  "glia" = c("SOX10", "S100B", "PLP1"),
  
  "Mast" = c("KIT", "TPSAB1", "TPSB2"),
  
  "MYLOID" = c("LYZ", "LST1", "S100A8"),
  
  "neuron" = c("SYT1", "SNAP25", "TUBB3"),
  
  "other" = c("MME", "CLU", "MMP7"),
  
  "pericyte" = c("ACTA2", "RGS5", "PDGFRB"),
  
  "plasma" = c("IGHA2", "JCHAIN", "MZB1"),
  
  "Stromal" = c("COL1A1", "DCN", "LUM"),
  
  "T/NK" = c("CD3E", "NKG7", "TRAC")
)

# =========================
# 3. 整理 marker 信息
# =========================
marker_df <- bind_rows(
  lapply(names(marker_genes_list), function(celltype_name) {
    data.frame(
      celltype = celltype_name,
      gene = marker_genes_list[[celltype_name]],
      stringsAsFactors = FALSE
    )
  })
)

marker_df <- marker_df %>%
  group_by(celltype) %>%
  mutate(gene_order = row_number()) %>%
  ungroup() %>%
  mutate(
    celltype = factor(celltype, levels = celltype_order)
  )

# =========================
# 4. 检查 marker 是否存在
# =========================
features_all <- unique(marker_df$gene)

features_use <- features_all[
  features_all %in% rownames(obj)
]

missing_genes <- setdiff(features_all, features_use)

if (length(missing_genes) > 0) {
  message("以下 marker 在 obj 中没有找到，已自动跳过：")
  message(paste(missing_genes, collapse = ", "))
}

marker_df <- marker_df %>%
  filter(gene %in% features_use)

# =========================
# 5. 计算 DotPlot 数据
# =========================
dp <- DotPlot(
  object = obj,
  features = features_use,
  group.by = "cell_subtype",
  scale = TRUE
)

df <- dp$data

# =========================
# 6. 合并 marker 对应关系
# =========================
df <- df %>%
  mutate(
    features.plot = as.character(features.plot),
    id = as.character(id)
  ) %>%
  inner_join(
    marker_df,
    by = c("features.plot" = "gene")
  )

# marker 按细胞类型顺序排列，但不分面
gene_order <- marker_df %>%
  arrange(celltype, gene_order) %>%
  pull(gene)

df$features.plot <- factor(
  df$features.plot,
  levels = gene_order
)

# y 轴细胞类型顺序
df$id <- factor(
  df$id,
  levels = rev(celltype_order)
)

# =========================
# 7. 配色
# =========================
orig_colors <- c(
  "#5749a0", "#0f7ab0", "#00bbb1", "#bef0b0", "#fdf4af",
  "#f9b64b", "#ec840e", "#ca443d", "#a51a49"
)

colors_14 <- colorRampPalette(orig_colors)(14)

breaks <- seq(0, 3, length.out = 15)

df$avg.exp.scaled.plot <- pmin(
  pmax(df$avg.exp.scaled, 0),
  3
)

# =========================
# 8. 不分层 DotPlot，一条线排开
# =========================
p <- ggplot(
  df,
  aes(
    x = features.plot,
    y = id,
    size = pct.exp,
    fill = avg.exp.scaled.plot
  )
) +
  geom_point(
    shape = 21,
    colour = "black",
    stroke = 0.5
  ) +
  scale_fill_stepsn(
    colours = colors_14,
    breaks = breaks,
    limits = c(0, 3),
    labels = NULL,
    show.limits = FALSE,
    name = "Avg. Exp. Scaled",
    guide = guide_colorsteps(
      title = "Avg. Exp. Scaled",
      barwidth = 1,
      barheight = 6,
      direction = "vertical",
      ticks = FALSE,
      frame.colour = "black",
      frame.linewidth = 1,
      label = FALSE
    )
  ) +
  scale_size_continuous(
    range = c(1, 8),
    limits = c(0, 100),
    name = "Percent Expressed"
  ) +
  theme(
    panel.background = element_blank(),
    panel.border = element_rect(
      fill = NA,
      color = "black",
      linewidth = 1.5
    ),
    panel.grid = element_blank(),
    
    axis.title = element_blank(),
    axis.text.y = element_text(
      color = "black",
      size = 12
    ),
    axis.text.x = element_text(
      color = "black",
      size = 12,
      angle = 90,
      hjust = 1,
      vjust = 0.5
    ),
    
    legend.position = "right",
    legend.box.background = element_blank(),
    legend.margin = margin(0, 0, 0, 10),
    legend.title = element_text(size = 10, color = "black"),
    legend.text = element_text(size = 9, color = "black")
  )

print(p)

# =========================
# 9. 保存 PDF
# =========================
ggsave(
  filename = "data/cell_subtype_marker_dotplot_single_line.pdf",
  plot = p,
  width = 16,
  height = 7,
  units = "in"
)

# =========================
# 9. 保存 PDF
# =========================
ggsave(
  filename = "data/cell_subtype_marker_dotplot_ordered.pdf",
  plot = p,
  width = 18,
  height = 8,
  units = "in"
)