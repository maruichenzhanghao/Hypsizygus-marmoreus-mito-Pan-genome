#!/usr/bin/env Rscript
# 内含子 presence/absence 热图 + 长度信息
# 行 = 内含子位点(按基因分组), 列 = 31个菌株(按系统发育聚类)

library(ggplot2)
library(ggsci)
library(reshape2)
library(cowplot)

setwd("/home/maxinxin/workspace/01.mito/03.pangenome.cactus/intron_analysis")

# ---- 读取数据 ----
pa <- read.delim("intron_pa_matrix.tsv", row.names = 1, check.names = FALSE)
len_mat <- read.delim("intron_length_matrix.tsv", row.names = 1, check.names = FALSE)
genome_stats <- read.delim("genome_intron_stats.tsv", stringsAsFactors = FALSE)

# ---- 准备内含子位点标签 ----
# 列名格式: cox1_P234, cox2_P324, etc.
col_labels <- colnames(pa)
gene_part <- gsub("_P.*", "", col_labels)
pos_part <- as.integer(gsub(".*_P", "", col_labels))

# 创建更友好的标签: cox1-i1(P234), cox1-i2(P380), ...
# 按基因内顺序编号
intron_labels <- c()
for (gene in unique(gene_part)) {
  idx <- which(gene_part == gene)
  for (j in seq_along(idx)) {
    intron_labels[idx[j]] <- paste0(gene, "-i", j, " (P", pos_part[idx[j]], ")")
  }
}

# ---- 菌株按系统发育关系聚类 ----
# 使用内含子模式聚类
dist_mat <- dist(pa, method = "binary")
hc <- hclust(dist_mat, method = "ward.D2")
strain_order <- rownames(pa)[hc$order]

# ---- 构建长表数据 ----
# 用长度值，0 = 缺失
len_long <- melt(as.matrix(len_mat))
colnames(len_long) <- c("strain", "intron", "length")
len_long$strain <- factor(len_long$strain, levels = strain_order)
len_long$intron <- factor(len_long$intron, levels = rev(col_labels))

# 添加基因分组
len_long$gene <- gsub("_P.*", "", as.character(len_long$intron))

# PA值
pa_long <- melt(as.matrix(pa))
colnames(pa_long) <- c("strain", "intron", "presence")
len_long$presence <- pa_long$presence

# 更友好的y轴标签
intron_label_map <- setNames(intron_labels, col_labels)
len_long$intron_label <- intron_label_map[as.character(len_long$intron)]
len_long$intron_label <- factor(len_long$intron_label, 
                                 levels = rev(intron_label_map))

# ---- 统计每个位点的频率 ----
freq_data <- data.frame(
  intron = col_labels,
  intron_label = intron_labels,
  freq = colSums(pa),
  gene = gene_part
)
freq_data$intron_label <- factor(freq_data$intron_label, levels = rev(intron_label_map))
freq_data$freq_label <- paste0(freq_data$freq, "/31")

# ---- 配色 ----
npg_colors <- pal_npg()(10)
gene_colors <- c("cox1" = npg_colors[1], "cox2" = npg_colors[2], 
                  "cob" = npg_colors[3], "nad5" = npg_colors[4])

# ---- 主热图 ----
p_main <- ggplot(len_long, aes(x = strain, y = intron_label)) +
  # 缺失的用浅灰色
  geom_tile(data = len_long[len_long$presence == 0, ], 
            fill = "grey95", color = "grey80", linewidth = 0.3) +
  # 存在的用蓝色渐变（按长度）
  geom_tile(data = len_long[len_long$presence == 1, ], 
            aes(fill = length), color = "grey40", linewidth = 0.3) +
  scale_fill_gradient(low = "#B3CDE3", high = "#023858",
                      name = "Intron\nlength (bp)",
                      limits = c(1000, 2600),
                      breaks = c(1000, 1500, 2000, 2500)) +
  # 基因分组色带 - 通过小面板分离
  labs(x = "", y = "") +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 10),
    axis.text.y = element_text(size = 11, face = "italic"),
    panel.grid = element_blank(),
    legend.position = "right",
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10),
    plot.margin = margin(5, 5, 5, 5)
  )

# ---- 右侧频率条形图 ----
p_freq <- ggplot(freq_data, aes(x = freq, y = intron_label, fill = gene)) +
  geom_col(width = 0.7, show.legend = TRUE) +
  geom_text(aes(label = freq_label), hjust = -0.1, size = 3.5) +
  scale_fill_manual(values = gene_colors, name = "Gene") +
  scale_x_continuous(limits = c(0, 38), breaks = c(0, 10, 20, 31)) +
  labs(x = "Frequency", y = "") +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "right",
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10, face = "italic"),
    plot.margin = margin(5, 5, 5, 0)
  )

# ---- 组合图 ----
p_combined <- plot_grid(p_main, p_freq, 
                        nrow = 1, rel_widths = c(3.5, 1.2),
                        align = "h", axis = "tb")

ggsave("plotD_intron_heatmap.pdf", p_combined, width = 16, height = 6, dpi = 300)
ggsave("plotD_intron_heatmap.png", p_combined, width = 16, height = 6, dpi = 300)
cat("plotD_intron_heatmap.pdf/png saved.\n")
