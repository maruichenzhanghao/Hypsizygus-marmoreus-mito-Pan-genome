# 加载必要的R包
setwd("/home/maruichen/workspace/9.R/02.mito/06.snpplote/")

library(ggplot2)   # 核心绘图包
library(dplyr)     # 数据处理
library(viridis)   # 科研级颜色方案

# 1. 输入数据（保持原始基因顺序）
df <- data.frame(
  gene = c("cox1", "nad4", "nad2", "nad3", "rps3", "nad6", "apt6", 
           "apt8", "nad5", "nad4L", "cox3", "cytb", "nad1", "cox2", "apt9"),
  pi = c(0.00382, 0.00128, 0.02533, 0.00317, 0.00217, 0.0019, 0.00129, 
         0, 0.00052, 0.0016, 0, 0.00410, 0.00012, 0.00059, 0.00145)
)

# 2. 锁定原始基因顺序
df <- df %>% 
  mutate(gene = factor(gene, levels = gene))

# 3. 绘制柱状图（添加外围框线）
p <- ggplot(df, aes(x = gene, y = pi, fill = pi)) +
  geom_col(width = 0.7, color = "black", linewidth = 0.3) +  # 柱状图主体
  
  # 颜色方案
  scale_fill_viridis(option = "plasma", name = "pi Value") +
  
  # 数值标签（根据数值大小调整位置）
  geom_text(
    aes(label = sprintf("%.4f", pi)),
    vjust = ifelse(df$pi == 0, 1.5, ifelse(df$pi < 0.005, 1.5, -0.5)),
    size = 3.2, 
    color = "#333333"
  ) +
  
  # 英文标题和坐标轴标签
  labs(
    title = NULL,
    x = "Gene", 
    y = "pi Value"
  ) +
  
  # 主题设置（添加外围框线）
  theme_minimal() +
  theme(
    # 添加外围框线
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.6),
    # 确保绘图区域背景透明，不与边框冲突
    panel.background = element_blank(),
    
    plot.title = element_text(
      hjust = 0.5, size = 14, face = "bold", margin = margin(b = 12)
    ),
    axis.title.x = element_text(size = 12, face = "bold", margin = margin(t = 10)),
    axis.title.y = element_text(size = 12, face = "bold", margin = margin(r = 10)),
    axis.text.x = element_text(
      angle = 45, hjust = 1, size = 10, color = "#555555"
    ),
    axis.text.y = element_text(size = 10, color = "#555555"),
    legend.title = element_text(size = 11, face = "bold"),
    legend.position = "right",
    plot.margin = margin(15, 15, 20, 15)
  ) +
  
  # 调整y轴范围
  ylim(0, max(df$pi) * 1.3)

# 显示图形
print(p)

# 保存为SVG矢量图
ggsave(
  "gene_pi_plot_with_border4.svg",
  plot = p,
  width = 9, height = 4,
  device = "svg"
)
