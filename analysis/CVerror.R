# 加载所需包
library(ggplot2)
library(ggsci)  # 科研级配色
library(extrafont)  # 可选，用于更丰富的字体

# 1. 准备数据
cv_data <- data.frame(
  K = 2:8,
  CV_error = c(0.72874, 0.50702, 0.73907, 0.62806, 0.51496, 0.12667, 0.21168)
)

# 确定最佳K值（CV误差最小的K）
best_k <- cv_data$K[which.min(cv_data$CV_error)]
best_error <- cv_data$CV_error[which.min(cv_data$CV_error)]

# 2. 绘制折线图
p <- ggplot(cv_data, aes(x = K, y = CV_error)) +
  # 折线设置：深蓝色粗线，带轻微阴影
  geom_line(color = "#0073C2FF", size = 1.2, lineend = "round") +
  # 数据点设置：白色填充，蓝色边框，带阴影
  geom_point(shape = 21, fill = "white", color = "#0073C2FF", 
             size = 4, stroke = 1.5, 
             position = position_dodge(0.1)) +
  # 标注最佳K值的垂直线
  geom_vline(xintercept = best_k, color = "#EFC000FF", 
             linetype = "dashed", size = 1) +
  # 标注最佳K值的文本
  annotate("text", x = best_k + 0.5, y = max(cv_data$CV_error) * 0.9, 
           label = paste0("Best K = ", best_k), 
           color = "#EFC000FF", fontface = "bold", size = 5) +
  # 坐标轴范围设置
  scale_x_continuous(limits = c(1.8, 8.2), breaks = 2:8, expand = c(0, 0)) +
  scale_y_continuous(limits = range(cv_data$CV_error) + c(-0.05, 0.05), expand = c(0, 0)) +
  # 标题和坐标轴标签
  labs(
    title = "Cross-Validation Error for Admixture Analysis",
    x = "Number of Populations (K)",
    y = "Cross-Validation Error"
  ) +
  # 主题美化
  theme_classic() +
  theme(
    # 标题设置
    plot.title = element_text(
      hjust = 0.5,  # 居中
      size = 16, 
      face = "bold", 
      color = "#333333",
      margin = margin(b = 15)
    ),
    # 坐标轴标签
    axis.title.x = element_text(
      size = 14, 
      face = "bold", 
      color = "#333333",
      margin = margin(t = 10)
    ),
    axis.title.y = element_text(
      size = 14, 
      face = "bold", 
      color = "#333333",
      margin = margin(r = 10)
    ),
    # 坐标轴刻度
    axis.text.x = element_text(size = 12, color = "#333333"),
    axis.text.y = element_text(size = 12, color = "#333333"),
    # 网格线（轻微显示，辅助读数）
    panel.grid.major.y = element_line(color = "#EEEEEE", size = 0.5),
    # 去除顶部和右侧边框
    axis.line = element_line(color = "#333333", size = 0.7)
  )

# 3. 保存图形（PDF格式适合期刊投稿）
ggsave(
  "cv_error_beautiful.pdf",
  plot = p,
  width = 8,
  height = 6,
  dpi = 300,
  device = "pdf"
)

# 也可以保存为PNG格式
ggsave(
  "cv_error_beautiful.png",
  plot = p,
  width = 8,
  height = 6,
  dpi = 300,
  device = "png"
)

# 显示图形
print(p)
