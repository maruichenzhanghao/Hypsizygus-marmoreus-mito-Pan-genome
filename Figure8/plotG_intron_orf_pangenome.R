#!/usr/bin/env Rscript
# 内含子ORF与泛基因组OG关联分析可视化
# Figure: cox1基因结构+ORF展示 + 内含子ORF OG分类汇总

library(ggplot2)
library(ggsci)
library(cowplot)
library(reshape2)

setwd("/home/maxinxin/workspace/01.mito/03.pangenome.cactus/intron_analysis")

# ---- 读取数据 ----
orf_data <- read.delim("intron_orf_og_full.tsv", stringsAsFactors = FALSE)
gs <- read.delim("gene_structure.tsv", stringsAsFactors = FALSE)

npg <- pal_npg()(10)

# ===== Panel A: cox1 基因结构 + 内嵌ORF =====
# 选择一个代表菌株展示完整的 cox1 外显子-内含子-ORF 结构
# Pattern 1 (24 strains): f2
rep_strain <- "f2"

cox1_struct <- gs[gs$strain == rep_strain & gs$gene == "cox1", ]
cox1_orfs <- orf_data[orf_data$strain == rep_strain & orf_data$host_gene == "cox1", ]

# 构建ORF在基因结构中的相对位置
# 内含子按顺序排列，ORF在每个内含子内
# 这里用简化版展示：在一条线上画出exon blocks + intron regions + ORF boxes

exons <- cox1_struct[cox1_struct$element_type == "exon", ]
introns <- cox1_struct[cox1_struct$element_type == "intron", ]

# 内含子ORF位置 (在每个内含子的中间位置)
intron_positions <- data.frame(
  intron_idx = introns$element_idx,
  intron_start = introns$rel_start,
  intron_end = introns$rel_end,
  intron_mid = (introns$rel_start + introns$rel_end) / 2,
  intron_len = introns$length
)

# 匹配ORF信息
# cox1 intron顺序对应的CDS位点
intron_cds_pos <- c(234, 380, 487, 609, 894, 1101, 1299)  # f2 pattern
orf_info <- data.frame(
  intron_idx = 1:7,
  cds_pos = intron_cds_pos,
  stringsAsFactors = FALSE
)

# 合并
intron_plot <- merge(intron_positions, orf_info, by = "intron_idx")

# 添加OG和类型信息
for (i in 1:nrow(intron_plot)) {
  pos <- intron_plot$cds_pos[i]
  sub <- cox1_orfs[cox1_orfs$intron_pos == pos, ]
  if (nrow(sub) > 0) {
    intron_plot$og[i] <- sub$og[1]
    intron_plot$orf_type[i] <- sub$orf_type[1]
    intron_plot$orf_name[i] <- sub$orf_name[1]
    intron_plot$og_n[i] <- sub$og_n_species[1]
    intron_plot$og_class[i] <- sub$og_class[1]
  }
}

# ORF type colors
orf_type_colors <- c("LAGLIDADG" = npg[1], "GIY-YIG" = npg[2], "Unknown" = "grey60")
og_class_colors <- c("Core" = npg[3], "Softcore" = npg[4], 
                      "Dispensable" = npg[5], "Unique" = npg[6])

max_x <- max(cox1_struct$rel_end)

p_a <- ggplot() +
  # 基因全长线
  geom_segment(aes(x = 0, xend = max_x, y = 0.5, yend = 0.5),
               color = "grey70", linewidth = 0.4) +
  # 外显子
  geom_rect(data = exons,
            aes(xmin = rel_start, xmax = rel_end, ymin = 0.25, ymax = 0.75),
            fill = npg[9], color = "grey30", linewidth = 0.4) +
  # 内含子区域 (浅色背景)
  geom_rect(data = intron_plot,
            aes(xmin = intron_start, xmax = intron_end, ymin = 0.35, ymax = 0.65),
            fill = "grey92", color = "grey70", linewidth = 0.3) +
  # ORF在内含子中 (彩色小方块)
  geom_rect(data = intron_plot,
            aes(xmin = intron_mid - intron_len * 0.3,
                xmax = intron_mid + intron_len * 0.3,
                ymin = 0.38, ymax = 0.62,
                fill = orf_type),
            color = "grey30", linewidth = 0.3, alpha = 0.85) +
  scale_fill_manual(values = orf_type_colors, name = "HEG type") +
  # 标注OG编号和CDS位点
  geom_text(data = intron_plot,
            aes(x = intron_mid, y = 0.85, 
                label = paste0(og, "\n(", og_class, ")")),
            size = 2.5, lineheight = 0.85) +
  geom_text(data = intron_plot,
            aes(x = intron_mid, y = 0.15,
                label = paste0("P", cds_pos)),
            size = 2.8, color = "grey40") +
  # 外显子编号
  geom_text(data = exons,
            aes(x = (rel_start + rel_end) / 2, y = 0.5,
                label = paste0("E", element_idx)),
            size = 2.5, color = "white", fontface = "bold") +
  labs(title = expression(italic("cox1") * " gene structure (f2, Pattern I)"),
       x = "Position (bp)", y = "") +
  scale_x_continuous(labels = function(x) paste0(round(x/1000, 1), "k")) +
  scale_y_continuous(limits = c(-0.05, 1.1)) +
  theme_minimal(base_size = 13) +
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    panel.grid = element_blank(),
    legend.position = c(0.92, 0.85),
    legend.background = element_rect(fill = alpha("white", 0.8)),
    plot.title = element_text(size = 14, hjust = 0.5)
  )

# ===== Panel B: 内含子ORF的OG分类在泛基因组中的角色 =====
# 每个OG一条，显示它在哪个内含子，OG频率，ORF类型

og_summary <- unique(orf_data[, c("og", "og_n_species", "og_class", "intron_label", "orf_type")])
# 合并同一个OG的多个内含子
og_agg <- aggregate(intron_label ~ og + og_n_species + og_class, 
                     data = og_summary, FUN = function(x) paste(unique(x), collapse = ", "))
og_agg_type <- aggregate(orf_type ~ og, data = og_summary, 
                          FUN = function(x) paste(unique(x), collapse = "/"))
og_agg <- merge(og_agg, og_agg_type, by = "og")
og_agg <- og_agg[order(-og_agg$og_n_species), ]

og_agg$og_label <- paste0(og_agg$og, " (", og_agg$intron_label, ")")
og_agg$og_label <- factor(og_agg$og_label, levels = rev(og_agg$og_label))
og_agg$og_class <- factor(og_agg$og_class, levels = c("Core", "Softcore", "Dispensable", "Unique"))

p_b <- ggplot(og_agg, aes(x = og_n_species, y = og_label, fill = og_class)) +
  geom_col(width = 0.7, color = "grey30", linewidth = 0.3) +
  geom_text(aes(label = orf_type), hjust = -0.1, size = 3.3, fontface = "italic") +
  scale_fill_manual(values = og_class_colors, name = "Pangenome\ncategory", drop = FALSE) +
  scale_x_continuous(limits = c(0, 40), breaks = c(0, 5, 10, 15, 20, 25, 30, 31)) +
  geom_vline(xintercept = c(2, 30, 31), linetype = "dashed", color = "grey60", linewidth = 0.3) +
  labs(x = "Number of strains with this OG", y = "",
       title = "Intron-encoded ORF orthogroups in pangenome") +
  theme_cowplot(font_size = 13) +
  theme(
    axis.text.y = element_text(size = 9),
    legend.position = c(0.75, 0.25),
    legend.background = element_rect(fill = alpha("white", 0.8)),
    plot.title = element_text(size = 14, hjust = 0.5)
  )

# ===== Panel C: 内含子ORF类型 vs 内含子频率 散点图 =====
# 每个内含子位点：频率 vs ORF长度，颜色=ORF类型

# 每个内含子位点一条记录
intron_summary <- aggregate(cbind(orf_len, og_n_species) ~ intron_label + orf_type + host_gene,
                            data = orf_data,
                            FUN = function(x) round(mean(x)))

# 获取每个内含子位点的出现频率（从PA矩阵）
pa <- read.delim("intron_pa_matrix.tsv", row.names = 1, check.names = FALSE)
intron_freq <- colSums(pa)

intron_summary$intron_freq <- NA
for (i in 1:nrow(intron_summary)) {
  label <- intron_summary$intron_label[i]
  # 转换: cox1-P234 -> cox1_P234
  pa_key <- gsub("-", "_", label)
  if (pa_key %in% names(intron_freq)) {
    intron_summary$intron_freq[i] <- intron_freq[pa_key]
  }
}

# 去重：同一个内含子+ORF类型只保留一条
intron_summary <- intron_summary[!is.na(intron_summary$intron_freq), ]

p_c <- ggplot(intron_summary, aes(x = intron_freq, y = orf_len, 
                                   color = orf_type, shape = host_gene)) +
  geom_point(size = 4, alpha = 0.85) +
  scale_color_manual(values = orf_type_colors, name = "HEG type") +
  scale_shape_manual(values = c("cox1" = 16, "cox2" = 17, "cob" = 15, "nad5" = 18),
                     name = "Host gene") +
  labs(x = "Intron frequency (n/31 strains)",
       y = "Mean ORF length (bp)",
       title = "Intron ORF length vs. intron frequency") +
  scale_x_continuous(breaks = c(1, 2, 5, 10, 15, 20, 25, 28, 30, 31)) +
  theme_cowplot(font_size = 13) +
  theme(
    legend.position = c(0.02, 0.98),
    legend.justification = c(0, 1),
    legend.background = element_rect(fill = alpha("white", 0.8)),
    plot.title = element_text(size = 14, hjust = 0.5)
  )

# ===== 组合 =====
p_top <- p_a
p_bottom <- plot_grid(p_b, p_c, nrow = 1, rel_widths = c(1.2, 1),
                      labels = c("B", "C"), label_size = 16)
p_final <- plot_grid(p_top, p_bottom, nrow = 2, 
                     rel_heights = c(0.7, 1.3),
                     labels = c("A", ""), label_size = 16)

ggsave("plotG_intron_orf_pangenome.pdf", p_final, width = 15, height = 11, dpi = 300)
ggsave("plotG_intron_orf_pangenome.png", p_final, width = 15, height = 11, dpi = 300)
cat("plotG_intron_orf_pangenome.pdf/png saved.\n")
