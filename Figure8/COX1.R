# 安装必要的包
if (!require("tidyverse")) install.packages("tidyverse")
if (!require("gggenes")) install.packages("gggenes")
if (!require("ggrepel")) install.packages("ggrepel")

library(tidyverse)
library(gggenes)
library(ggrepel)

# 读取数据
gene_data <- read.csv("cox1_features.xiugai.csv")

# 查看数据结构
head(gene_data)
str(gene_data)

# 创建颜色映射
type_colors <- c(
  "gene" = "#2c3e50",
  "exon" = "#e74c3c",
  "intron" = "#3498db",
  "orf" = "#27ae60"
)

# 绘制单个样本的基因结构图
plot_individual_gene <- function(sample_name, data) {
  sample_data <- data %>% filter(molecule == sample_name)
  
  p <- ggplot(sample_data, 
              aes(xmin = start, xmax = end, y = molecule, 
                  fill = type, label = gene)) +
    # 绘制基因区域
    geom_gene_arrow(arrowhead_height = unit(3, "mm"), 
                    arrowhead_width = unit(1, "mm")) +
    # 添加标签
    geom_gene_label(align = "left", grow = TRUE) +
    # 设置颜色
    scale_fill_manual(values = type_colors) +
    # 主题设置
    theme_genes() +
    labs(title = paste("cox1 Gene Structure -", sample_name),
         x = "Position (bp)",
         y = "",
         fill = "Region Type") +
    theme(legend.position = "bottom",
          plot.title = element_text(hjust = 0.5, face = "bold"),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank())
  
  return(p)
}

# 为每个样本绘制单独的图
samples <- unique(gene_data$molecule)

# 绘制并保存每个样本的图
for (sample in samples) {
  p <- plot_individual_gene(sample, gene_data)
  ggsave(paste0("gene_structure_", sample, ".png"), p, 
         width = 12, height = 6, dpi = 300)
  print(p)  # 在R中显示图
}

# 绘制所有样本的比较图（简化版）
plot_comparison_simple <- function(data) {
  # 只显示主要区域类型
  main_data <- data %>% filter(type %in% c("exon", "intron"))
  
  p <- ggplot(main_data, 
              aes(xmin = start, xmax = end, y = molecule, 
                  fill = type, forward = TRUE)) +
    geom_gene_arrow() +
    scale_fill_manual(values = c("exon" = "#e74c3c", "intron" = "#3498db")) +
    theme_genes() +
    labs(title = "Comparative cox1 Gene Structure",
         x = "Position (bp)",
         y = "Sample",
         fill = "Region Type") +
    theme(legend.position = "bottom",
          axis.text.y = element_text(size = 8),
          plot.title = element_text(hjust = 0.5, face = "bold"))
  
  ggsave("cox1_comparison_simple.png", p, width = 14, height = 10, dpi = 300)
  return(p)
}

# 绘制详细的比较图（分组显示）
plot_comparison_detailed <- function(data) {
  p <- ggplot(data, 
              aes(xmin = start, xmax = end, y = molecule, 
                  fill = type, forward = TRUE)) +
    geom_gene_arrow() +
    facet_grid(type ~ ., scales = "free_y", space = "free") +
    scale_fill_manual(values = type_colors) +
    theme_genes() +
    labs(title = "Detailed cox1 Gene Structure Comparison",
         x = "Position (bp)",
         y = "Sample",
         fill = "Region Type") +
    theme(legend.position = "none",
          strip.text.y = element_text(angle = 0),
          axis.text.y = element_text(size = 6))
  
  ggsave("cox1_comparison_detailed.png", p, width = 16, height = 12, dpi = 300)
  return(p)
}

# 绘制ORF区域的特殊图
plot_orf_regions <- function(data) {
  orf_data <- data %>% filter(type == "orf")
  
  p <- ggplot(orf_data, 
              aes(xmin = start, xmax = end, y = molecule, 
                  fill = gene, label = gene)) +
    geom_gene_arrow() +
    geom_gene_label(align = "left") +
    scale_fill_brewer(palette = "Set3") +
    theme_genes() +
    labs(title = "ORF Regions within cox1 Gene",
         x = "Position (bp)",
         y = "Sample",
         fill = "ORF Type") +
    theme(legend.position = "bottom",
          plot.title = element_text(hjust = 0.5, face = "bold"))
  
  ggsave("orf_regions.png", p, width = 14, height = 8, dpi = 300)
  return(p)
}

# 执行绘图
comparison_simple <- plot_comparison_simple(gene_data)
comparison_detailed <- plot_comparison_detailed(gene_data)
orf_plot <- plot_orf_regions(gene_data)

# 显示图形
print(comparison_simple)
print(comparison_detailed)
print(orf_plot)






# 高级绘图：带缩放和更好的标签
plot_advanced_gene_structure <- function(data) {
  # 计算每个样本的总长度
  sample_lengths <- data %>%
    group_by(molecule) %>%
    summarize(total_length = max(end))
  
  # 合并数据
  plot_data <- data %>%
    left_join(sample_lengths, by = "molecule") %>%
    mutate(relative_start = start / total_length,
           relative_end = end / total_length)
  
  p <- ggplot(plot_data, 
              aes(xmin = relative_start, xmax = relative_end, y = molecule, 
                  fill = type, alpha = type)) +
    geom_gene_arrow() +
    scale_fill_manual(values = type_colors) +
    scale_alpha_manual(values = c("gene" = 0.3, "exon" = 1.0, 
                                  "intron" = 0.6, "orf" = 0.8)) +
    theme_genes() +
    labs(title = "Advanced cox1 Gene Structure Comparison (Scaled)",
         x = "Relative Position",
         y = "Sample",
         fill = "Region Type",
         alpha = "Region Type") +
    theme(legend.position = "bottom",
          plot.title = element_text(hjust = 0.5, face = "bold"))
  
  ggsave("cox1_advanced_comparison.png", p, width = 14, height = 10, dpi = 300)
  return(p)
}

# 执行高级绘图
advanced_plot <- plot_advanced_gene_structure(gene_data)
print(advanced_plot)
