library(ggplot2)

# 生成示例数据
set.seed(42)
df <- data.frame(
  features.plot = rep(c("Feature1", "Feature2", "Feature3", "Feature4", "Feature5"), times = 5),
  id = rep(1:5, each = 5),
  celltype = rep(c("CellTypeA", "CellTypeB", "CellTypeC", "CellTypeD", "CellTypeE"), times = 5),
  pct.exp = runif(25, 5, 100),
  avg.exp.scaled = runif(25, 0, 3)
)

# 定义原始9色，扩展为14色
orig_colors <- c('#5749a0', '#0f7ab0', '#00bbb1', '#bef0b0', '#fdf4af', 
                 '#f9b64b', '#ec840e', '#ca443d', '#a51a49')
colors_14 <- colorRampPalette(orig_colors)(14)

# 设定14个区间的断点（0-3等分为15个点）
breaks <- seq(0, 3, length.out = 15)

# 绘图
p <- ggplot(df, aes(features.plot, interaction(id, celltype),
                    size = pct.exp, fill = avg.exp.scaled)) + 
  geom_point(shape = 21, colour = "black", stroke = 0.5) +
  theme(
    panel.background = element_blank(),
    # 加粗面板边框（主图框架）
    panel.border = element_rect(fill = NA, color = "black", linewidth = 1.5),
    panel.grid = element_blank(),
    axis.title = element_blank(),
    axis.text.y = element_text(color = 'black', size = 12),
    axis.text.x = element_text(color = 'black', size = 12, angle = 90, hjust = 1, vjust = 0.5),
    legend.position = "right",
    legend.box.background = element_blank(),
    legend.margin = margin(0, 0, 0, 10)
  ) +
  scale_fill_stepsn(
    colours = colors_14,
    breaks = breaks,
    # 核心修改1：移除所有颜色标签
    labels = NULL,
    # 核心修改2：关闭首尾标签显示
    show.limits = FALSE,
    guide = guide_colorsteps(
      title = "Avg. Exp. Scaled",
      barwidth = 1,
      barheight = 6,
      direction = "vertical",
      ticks = FALSE,
      frame.colour = "black",
      frame.linewidth = 1,
      # 核心修改3：确保颜色条不显示任何标签
      label = FALSE
    )
  ) +
  # 补充size的刻度范围（可选）
  scale_size_continuous(range = c(1, 8))

print(p)