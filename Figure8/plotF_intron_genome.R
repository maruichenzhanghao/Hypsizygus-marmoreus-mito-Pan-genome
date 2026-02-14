#!/usr/bin/env Rscript
# 内含子数量与基因组大小的相关性分析
# 散点图 + 回归 + 内含子/外显子组成堆叠条形图

library(ggplot2)
library(ggsci)
library(cowplot)

setwd("/home/maxinxin/workspace/01.mito/03.pangenome.cactus/intron_analysis")

# ---- 读取数据 ----
stats <- read.delim("genome_intron_stats.tsv", stringsAsFactors = FALSE)
summary_data <- read.delim("strain_intron_summary.tsv", stringsAsFactors = FALSE)

# 合并 - 避免列名冲突
summary_data2 <- summary_data[, c("strain", "cox1_introns", "cox2_introns", "cob_introns", "nad5_introns")]
df <- merge(stats, summary_data2, by = "strain")

# 确保数值列为数值类型
df$genome_len <- as.numeric(df$genome_len)
df$total_introns <- as.numeric(df$total_introns)
df$total_intron_len <- as.numeric(df$total_intron_len)

# 按基因组大小排序
df <- df[order(df$genome_len), ]

# 定义菌株分组
df$group <- "Majority (n=24)"
df$group[df$strain %in% c("SRR12151860", "SRR12151871", "SRR12151883")] <- "Rearranged (n=3)"
df$group[df$strain %in% c("nn12-1", "nn12-17")] <- "nn12 (n=2)"
df$group[df$strain == "SRR7874787"] <- "SRR7874787"
df$group[df$strain == "MH382825.1"] <- "MH382825.1"

npg_colors <- pal_npg()(10)
group_colors <- c(
  "Majority (n=24)" = npg_colors[1],
  "Rearranged (n=3)" = npg_colors[2],
  "nn12 (n=2)" = npg_colors[3],
  "SRR7874787" = npg_colors[4],
  "MH382825.1" = npg_colors[5]
)

group_shapes <- c(
  "Majority (n=24)" = 16,
  "Rearranged (n=3)" = 17,
  "nn12 (n=2)" = 15,
  "SRR7874787" = 18,
  "MH382825.1" = 8
)

# ---- 图A: 散点图 - 内含子数 vs 基因组大小 ----
cor_test <- cor.test(df$total_introns, df$genome_len)
r_val <- round(cor_test$estimate, 3)
p_val <- format(cor_test$p.value, digits = 3, scientific = TRUE)

p_scatter <- ggplot(df, aes(x = total_introns, y = genome_len / 1000, 
                             color = group, shape = group)) +
  geom_point(size = 4, alpha = 0.85) +
  geom_smooth(aes(group = 1), method = "lm", se = TRUE, 
              color = "grey50", fill = "grey85", linewidth = 0.8, linetype = "dashed") +
  scale_color_manual(values = group_colors, name = "Group") +
  scale_shape_manual(values = group_shapes, name = "Group") +
  annotate("text", x = 7.5, y = 112, 
           label = paste0("R = ", r_val, "\np = ", p_val),
           size = 4.5, hjust = 0, fontface = "italic") +
  labs(x = "Number of introns", y = "Genome size (kb)",
       title = "Intron number vs. mitogenome size") +
  scale_x_continuous(breaks = 7:10) +
  theme_cowplot(font_size = 14) +
  theme(
    legend.position = c(0.02, 0.98),
    legend.justification = c(0, 1),
    legend.background = element_rect(fill = alpha("white", 0.8)),
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10),
    plot.title = element_text(size = 15, hjust = 0.5)
  )

# ---- 图B: 基因组组成堆叠条形图(按基因组大小排序) ----
df$coding_len <- df$genome_len - df$total_intron_len
# 进一步分解：外显子+其他编码 vs 内含子 vs 基因间区
# 简化为: 内含子 vs 非内含子
df_comp <- data.frame(
  strain = rep(df$strain, 2),
  component = rep(c("Intron", "Non-intron"), each = nrow(df)),
  length = c(df$total_intron_len, df$genome_len - df$total_intron_len),
  genome_len = rep(df$genome_len, 2)
)
df_comp$strain <- factor(df_comp$strain, levels = df$strain[order(df$genome_len)])
df_comp$component <- factor(df_comp$component, levels = c("Non-intron", "Intron"))

p_bar <- ggplot(df_comp, aes(x = strain, y = length / 1000, fill = component)) +
  geom_col(width = 0.75) +
  scale_fill_manual(values = c("Non-intron" = npg_colors[9], "Intron" = npg_colors[1]),
                    name = "Component") +
  labs(x = "", y = "Length (kb)", 
       title = "Mitogenome composition") +
  theme_cowplot(font_size = 14) +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 8),
    legend.position = c(0.02, 0.98),
    legend.justification = c(0, 1),
    legend.background = element_rect(fill = alpha("white", 0.8)),
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10),
    plot.title = element_text(size = 15, hjust = 0.5)
  )

# ---- 图C: 各基因内含子数量的lollipop图 ----
# 每个菌株每个基因的内含子数量
intron_by_gene <- data.frame(
  strain = rep(summary_data$strain, 4),
  gene = rep(c("cox1", "cox2", "cob", "nad5"), each = nrow(summary_data)),
  n_introns = c(summary_data$cox1_introns, summary_data$cox2_introns,
                summary_data$cob_introns, summary_data$nad5_introns)
)

# 箱线图汇总
gene_colors_v <- c("cox1" = npg_colors[1], "cox2" = npg_colors[2],
                    "cob" = npg_colors[3], "nad5" = npg_colors[4])

p_box <- ggplot(intron_by_gene, aes(x = gene, y = n_introns, fill = gene)) +
  geom_boxplot(width = 0.5, alpha = 0.7, outlier.shape = NA) +
  geom_jitter(width = 0.15, alpha = 0.5, size = 1.5, color = "grey30") +
  scale_fill_manual(values = gene_colors_v) +
  labs(x = "", y = "Number of introns",
       title = "Intron count per gene") +
  theme_cowplot(font_size = 14) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(face = "italic", size = 13),
    plot.title = element_text(size = 15, hjust = 0.5)
  )

# ---- 组合 ----
p_top_row <- plot_grid(p_scatter, p_box, nrow = 1, rel_widths = c(1.3, 1),
                       labels = c("A", "B"), label_size = 16)
p_final <- plot_grid(p_top_row, p_bar, nrow = 2, rel_heights = c(1, 0.9),
                     labels = c("", "C"), label_size = 16)

ggsave("plotF_intron_genome.pdf", p_final, width = 14, height = 10, dpi = 300)
ggsave("plotF_intron_genome.png", p_final, width = 14, height = 10, dpi = 300)
cat("plotF_intron_genome.pdf/png saved.\n")
