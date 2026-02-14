#!/usr/bin/env Rscript
# 基因结构对比图 - 展示4个基因的外显子-内含子结构变异
# 每个基因的不同模式并排展示，外显子用彩色方块，内含子用细线/浅色方块

library(ggplot2)
library(ggsci)
library(cowplot)

setwd("/home/maxinxin/workspace/01.mito/03.pangenome.cactus/intron_analysis")

# ---- 读取基因结构数据 ----
gs <- read.delim("gene_structure.tsv", stringsAsFactors = FALSE)

# ---- 定义每个基因的模式及其代表菌株 ----
# cox1: 5 patterns
# P1 (24): (234,380,487,609,894,1101,1299) - 7 introns
# P2 (3): (234,380,609,700,861,1101,1299) - SRR12151860/71/83
# P3 (2): (234,380,700,861,1101,1299) - nn12-1/17
# P4 (1): (234,380,609,894,1101,1299) - MH382825.1
# P5 (1): (234,380,703,861,1101,1299) - SRR7874787

# cox2: 2 patterns
# P1 (30): (324) - 1 intron
# P2 (1): () - SRR7874787

# cob: 3 patterns
# P1 (28): (490) - 1 intron
# P2 (2): (490,820) - nn12-1/17
# P3 (1): () - SRR7874787

# nad5: 2 patterns
# P1 (28): (324) - 1 intron
# P2 (3): () - SRR12151860/71/83

patterns <- list(
  cox1 = list(
    list(name = "Pattern I (n=24)", rep = "f2"),
    list(name = "Pattern II (n=3)", rep = "SRR12151860"),
    list(name = "Pattern III (n=2)", rep = "nn12-1"),
    list(name = "Pattern IV (n=1)", rep = "MH382825.1"),
    list(name = "Pattern V (n=1)", rep = "SRR7874787")
  ),
  cox2 = list(
    list(name = "Pattern I (n=30)", rep = "f2"),
    list(name = "Pattern II (n=1)", rep = "SRR7874787")
  ),
  cob = list(
    list(name = "Pattern I (n=28)", rep = "f2"),
    list(name = "Pattern II (n=2)", rep = "nn12-1"),
    list(name = "Pattern III (n=1)", rep = "SRR7874787")
  ),
  nad5 = list(
    list(name = "Pattern I (n=28)", rep = "f2"),
    list(name = "Pattern II (n=3)", rep = "SRR12151860")
  )
)

npg_colors <- pal_npg()(10)

# ---- 构建绘图数据 ----
plot_data <- data.frame()
label_data <- data.frame()

y_pos <- 0  # 从上到下递增
gene_label_y <- c()

for (gene in c("cox1", "cox2", "cob", "nad5")) {
  gene_patterns <- patterns[[gene]]
  gene_start_y <- y_pos
  
  for (p in gene_patterns) {
    y_pos <- y_pos + 1
    pname <- p$name
    rep_strain <- p$rep
    
    # 获取该菌株该基因的结构数据
    sub <- gs[gs$strain == rep_strain & gs$gene == gene, ]
    
    if (nrow(sub) == 0) next
    
    for (i in 1:nrow(sub)) {
      row <- sub[i, ]
      plot_data <- rbind(plot_data, data.frame(
        gene = gene,
        pattern = pname,
        y = y_pos,
        xmin = row$rel_start,
        xmax = row$rel_end,
        type = row$element_type,
        length = row$length,
        idx = row$element_idx,
        stringsAsFactors = FALSE
      ))
    }
    
    label_data <- rbind(label_data, data.frame(
      gene = gene,
      pattern = pname,
      y = y_pos,
      stringsAsFactors = FALSE
    ))
  }
  
  gene_label_y[[gene]] <- (gene_start_y + y_pos) / 2 + 0.5
  y_pos <- y_pos + 0.5  # gap between genes
}

# ---- 绘制每个基因的结构图 ----
make_gene_plot <- function(gene_name, gene_color) {
  sub_plot <- plot_data[plot_data$gene == gene_name, ]
  sub_label <- label_data[label_data$gene == gene_name, ]
  
  if (nrow(sub_plot) == 0) return(NULL)
  
  # 外显子和内含子分开
  exons <- sub_plot[sub_plot$type == "exon", ]
  introns <- sub_plot[sub_plot$type == "intron", ]
  
  # 计算最大范围
  max_x <- max(sub_plot$xmax)
  
  p <- ggplot() +
    # 中间连线 (基因全长线)
    geom_segment(data = sub_label, 
                 aes(x = 0, xend = max_x, y = y, yend = y),
                 color = "grey60", linewidth = 0.5) +
    # 内含子 - 用浅色填充的薄矩形
    geom_rect(data = introns,
              aes(xmin = xmin, xmax = xmax, ymin = y - 0.2, ymax = y + 0.2),
              fill = "grey85", color = "grey60", linewidth = 0.3) +
    # 外显子 - 用深色填充的厚矩形
    geom_rect(data = exons,
              aes(xmin = xmin, xmax = xmax, ymin = y - 0.35, ymax = y + 0.35),
              fill = gene_color, color = "grey30", linewidth = 0.4, alpha = 0.9) +
    # 内含子长度标签
    geom_text(data = introns,
              aes(x = (xmin + xmax) / 2, y = y + 0.38, 
                  label = paste0(round(length/1000, 1), "k")),
              size = 2.8, color = "grey40") +
    # Y轴标签
    scale_y_continuous(breaks = sub_label$y, labels = sub_label$pattern,
                       expand = c(0.15, 0.15)) +
    scale_x_continuous(labels = function(x) paste0(x/1000, "k"),
                       expand = c(0.02, 0.02)) +
    labs(title = gene_name, x = "Position (bp)", y = "") +
    theme_minimal(base_size = 13) +
    theme(
      plot.title = element_text(face = "italic", size = 15, hjust = 0.5),
      axis.text.y = element_text(size = 10),
      axis.text.x = element_text(size = 9),
      axis.title.x = element_text(size = 11),
      panel.grid.major.y = element_blank(),
      panel.grid.minor = element_blank(),
      plot.margin = margin(8, 10, 5, 5)
    )
  
  return(p)
}

p_cox1 <- make_gene_plot("cox1", npg_colors[1])
p_cox2 <- make_gene_plot("cox2", npg_colors[2])
p_cob <- make_gene_plot("cob", npg_colors[3])
p_nad5 <- make_gene_plot("nad5", npg_colors[4])

# ---- 组合图 ----
# cox1最复杂（5个pattern），其他较简单
p_top <- p_cox1
p_bottom <- plot_grid(p_nad5, p_cob, p_cox2, 
                      nrow = 1, rel_widths = c(1, 1, 1),
                      labels = c("", "", ""), align = "h")

p_final <- plot_grid(p_top, p_bottom, 
                     nrow = 2, rel_heights = c(1.5, 1),
                     labels = c("", ""))

ggsave("plotE_gene_structure.pdf", p_final, width = 14, height = 9, dpi = 300)
ggsave("plotE_gene_structure.png", p_final, width = 14, height = 9, dpi = 300)
cat("plotE_gene_structure.pdf/png saved.\n")
