library(ggplot2)
library(ggsci)
library(dplyr)
library(readr)
setwd("/home/maruichen/workspace/9.R/02.mito/03.kaks/")
# 读取数据
data <- read_csv("cleaned1.csv")

# 优化后的箱线图
p <- ggplot(data, aes(x = gene, y = omega, fill = gene)) +
  geom_boxplot(
    width = 0.6,
    alpha = 0.8,
    outlier.shape = 21,
    outlier.size = 2.5,
    outlier.color = "black",
    outlier.fill = "#E64B35"
  ) +
  geom_jitter(
    width = 0.2,
    size = 1.5,
    alpha = 0.5,
    color = "gray40"
  ) +
  scale_fill_npg() +
  labs(
    title = "Ka/Ks (ω) Ratio Distribution by Gene",
    x = "Gene",
    y = "ω (Ka/Ks)",
    caption = "Horizontal line indicates neutral evolution (ω=1)"
  ) +
  geom_hline(
    yintercept = 1,
    linetype = "dashed",
    color = "red",
    linewidth = 0.8
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      vjust = 1,
      face = "italic"
    ),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "none",
    plot.caption = element_text(hjust = 0.5, color = "gray40"),
    # 添加外围框线关键设置
    panel.border = element_rect(colour = "black", fill = NA, size = 0.8)
  ) +
  stat_summary(
    fun = median,
    geom = "text",
    aes(label = sprintf("%.3f", round(..y.., 3))),
    vjust = -0.8,
    size = 3.5,
    fontface = "bold",
    color = "black"
  ) +
  scale_y_continuous(breaks = seq(0, max(data$omega)+0.5, by = 0.5))

print(p)
# 保存为高质量PDF（矢量图）
ggsave(
  "omega_boxplot1.pdf", 
  plot = p,
  width = 12,             # 加宽以适应基因名
  height = 7,
  device = cairo_pdf      # 确保字体嵌入
)

# 可选：同时保存PNG版本
ggsave(
  "omega_boxplot.png",
  plot = p,
  width = 12,
  height = 7,
  dpi = 600,
  bg = "white"
)
