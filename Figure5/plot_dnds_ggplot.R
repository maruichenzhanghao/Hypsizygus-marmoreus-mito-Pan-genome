#!/usr/bin/env Rscript
# dN/dS分析结果可视化 - 使用ggplot2 + ggsci
# 真姬菇线粒体基因组15个核心基因

library(ggplot2)
library(ggsci)
library(dplyr)
library(tidyr)
library(cowplot)

# 读取数据
df <- read.csv("dnds_data_for_R.csv", stringsAsFactors = FALSE)

# 转换has_data为逻辑值
df$has_data <- df$has_data == "True"

# 基因顺序（按功能分组）
gene_order <- c('atp6', 'atp8', 'atp9', 'cox1', 'cox2', 'cox3', 'cob', 
                'nad1', 'nad2', 'nad3', 'nad4', 'nad4l', 'nad5', 'nad6', 'rps3')

# 只保留有数据的基因
df_valid <- df %>% filter(has_data == TRUE & !is.na(omega))
df_valid$gene <- factor(df_valid$gene, levels = gene_order)

cat("有效数据行数:", nrow(df_valid), "\n")
cat("包含的基因:", unique(as.character(df_valid$gene)), "\n\n")

# 计算每个基因的统计量
stats <- df_valid %>%
  group_by(gene) %>%
  summarise(
    n = n(),
    mean = mean(omega, na.rm = TRUE),
    median = median(omega, na.rm = TRUE),
    sd = sd(omega, na.rm = TRUE),
    .groups = 'drop'
  )

cat("各基因统计量:\n")
print(as.data.frame(stats))

# 定义配色
gene_colors <- c(
  'cox1' = '#4DBBD5', 'cob' = '#00A087', 'nad2' = '#3C5488',
  'nad3' = '#7E6148', 'nad4l' = '#8491B4', 'nad5' = '#91D1C2', 'rps3' = '#F39B7F'
)

# ============================================
# 图1: 小提琴图 + 箱线图 + 散点
# ============================================
p1 <- ggplot(df_valid, aes(x = gene, y = omega, fill = gene)) +
  geom_violin(alpha = 0.6, width = 0.8, trim = FALSE) +
  geom_boxplot(width = 0.15, fill = "white", alpha = 0.8, outlier.shape = NA) +
  geom_jitter(width = 0.1, size = 1, alpha = 0.4, color = "black") +
  geom_hline(yintercept = 1, linetype = "dashed", color = "red", linewidth = 1) +
  scale_fill_npg() +
  scale_y_continuous(limits = c(0, 1.1), breaks = seq(0, 1, 0.2)) +
  labs(
    title = "dN/dS (ω) Distribution by Gene",
    subtitle = "Hypsizygus marmoreus Mitochondrial Core Genes (Non-zero ω values)",
    x = "Gene",
    y = "ω (dN/dS)",
    caption = "Red dashed line: Neutral evolution (ω=1)\nOnly genes with valid dS>0 are shown"
  ) +
  theme_bw(base_size = 12) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
    plot.subtitle = element_text(hjust = 0.5, size = 11, color = "gray40"),
    axis.text.x = element_text(angle = 45, hjust = 1, face = "italic", size = 10),
    legend.position = "none",
    panel.grid.minor = element_blank()
  )

# ============================================
# 图2: 雨云图风格（半小提琴）
# ============================================
p2 <- ggplot(df_valid, aes(x = gene, y = omega, fill = gene)) +
  geom_violin(alpha = 0.7, width = 0.9) +
  geom_boxplot(width = 0.12, alpha = 0.9, outlier.size = 1.5, outlier.alpha = 0.6) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "#DC0000", linewidth = 0.8) +
  scale_fill_lancet() +
  scale_y_continuous(limits = c(0, 1.1), breaks = seq(0, 1, 0.2)) +
  labs(
    title = "Selection Pressure on Mitochondrial Genes",
    subtitle = "ω < 1: Purifying selection | ω = 1: Neutral | ω > 1: Positive selection",
    x = "Gene",
    y = "ω (dN/dS)"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
    plot.subtitle = element_text(hjust = 0.5, size = 10, color = "gray50"),
    axis.text.x = element_text(angle = 45, hjust = 1, face = "italic"),
    legend.position = "none",
    panel.grid.major.x = element_blank()
  )

# ============================================
# 图3: 条形图 + 误差棒（展示中位数±SD）
# ============================================
stats$sd[is.na(stats$sd)] <- 0

p3 <- ggplot(stats, aes(x = reorder(gene, -median), y = median, fill = gene)) +
  geom_bar(stat = "identity", alpha = 0.85, width = 0.7) +
  geom_errorbar(aes(ymin = pmax(0, median - sd), ymax = pmin(1.1, median + sd)), 
                width = 0.2, color = "gray30") +
  geom_point(size = 3, color = "black") +
  geom_hline(yintercept = 1, linetype = "dashed", color = "red", linewidth = 0.8) +
  geom_text(aes(label = sprintf("%.2f", median)), vjust = -0.5, size = 3.5) +
  scale_fill_jco() +
  scale_y_continuous(limits = c(0, 1.2), breaks = seq(0, 1, 0.2)) +
  labs(
    title = "Median ω Values with Standard Deviation",
    subtitle = "Genes ordered by ω value (descending)",
    x = "Gene",
    y = "Median ω (dN/dS)"
  ) +
  theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5, color = "gray50"),
    axis.text.x = element_text(angle = 45, hjust = 1, face = "italic"),
    legend.position = "none"
  )

# ============================================
# 图4: 15个基因综合展示（标注无数据的基因）
# ============================================
all_genes_stats <- data.frame(gene = gene_order, stringsAsFactors = FALSE)
all_genes_stats <- left_join(all_genes_stats, stats, by = "gene")
all_genes_stats$median[is.na(all_genes_stats$median)] <- 0
all_genes_stats$n[is.na(all_genes_stats$n)] <- 0
all_genes_stats$has_data <- all_genes_stats$n > 0

all_genes_stats$label <- ifelse(
  all_genes_stats$has_data,
  sprintf("%.2f\n(n=%d)", all_genes_stats$median, all_genes_stats$n),
  "N/A\n(dS=0)"
)

all_genes_stats$gene <- factor(all_genes_stats$gene, levels = gene_order)

# 功能分组
all_genes_stats$group <- case_when(
  all_genes_stats$gene %in% c('atp6', 'atp8', 'atp9') ~ "ATP synthase",
  all_genes_stats$gene %in% c('cox1', 'cox2', 'cox3') ~ "Cytochrome c oxidase",
  all_genes_stats$gene == 'cob' ~ "Cytochrome b",
  all_genes_stats$gene %in% c('nad1', 'nad2', 'nad3', 'nad4', 'nad4l', 'nad5', 'nad6') ~ "NADH dehydrogenase",
  all_genes_stats$gene == 'rps3' ~ "Ribosomal protein",
  TRUE ~ "Other"
)

p4 <- ggplot(all_genes_stats, aes(x = gene, y = median, fill = group)) +
  geom_bar(stat = "identity", alpha = 0.85, width = 0.75) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "red", linewidth = 0.8) +
  geom_text(aes(label = label, y = median + 0.08), size = 2.6, vjust = 0) +
  scale_fill_jco() +
  scale_y_continuous(limits = c(0, 1.3), breaks = seq(0, 1, 0.2)) +
  labs(
    title = "Selection Pressure on 15 Core Mitochondrial Genes",
    subtitle = "Hypsizygus marmoreus (31 strains)",
    x = "Gene",
    y = "Median ω (dN/dS)",
    fill = "Gene Function",
    caption = "N/A: Insufficient synonymous substitutions (dS≈0) to calculate ω"
  ) +
  theme_bw(base_size = 11) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 13),
    plot.subtitle = element_text(hjust = 0.5, color = "gray40"),
    axis.text.x = element_text(angle = 45, hjust = 1, face = "italic", size = 9),
    legend.position = "bottom",
    legend.title = element_text(face = "bold"),
    panel.grid.minor = element_blank()
  ) +
  guides(fill = guide_legend(nrow = 2))

# 保存图表
ggsave("dnds_violin_plot.png", p1, width = 10, height = 6, dpi = 300)
ggsave("dnds_violin_plot.pdf", p1, width = 10, height = 6)

ggsave("dnds_lancet_plot.png", p2, width = 10, height = 6, dpi = 300)
ggsave("dnds_lancet_plot.pdf", p2, width = 10, height = 6)

ggsave("dnds_barplot.png", p3, width = 9, height = 6, dpi = 300)
ggsave("dnds_barplot.pdf", p3, width = 9, height = 6)

ggsave("dnds_all_genes.png", p4, width = 11, height = 7, dpi = 300)
ggsave("dnds_all_genes.pdf", p4, width = 11, height = 7)

# 组合图
combined <- plot_grid(p1, p4, ncol = 1, labels = c("A", "B"), label_size = 14)
ggsave("dnds_combined.png", combined, width = 11, height = 12, dpi = 300)
ggsave("dnds_combined.pdf", combined, width = 11, height = 12)

cat("\n========================================\n")
cat("图表已成功保存:\n")
cat("  - dnds_violin_plot.png/pdf (小提琴图+箱线图)\n")
cat("  - dnds_lancet_plot.png/pdf (Lancet配色小提琴图)\n")
cat("  - dnds_barplot.png/pdf (条形图+误差棒)\n")
cat("  - dnds_all_genes.png/pdf (15基因综合图)\n")
cat("  - dnds_combined.png/pdf (组合图)\n")
cat("========================================\n")
