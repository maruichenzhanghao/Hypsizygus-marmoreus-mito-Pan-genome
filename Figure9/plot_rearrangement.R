#!/usr/bin/env Rscript
# =============================================================================
# 白玉菇线粒体基因组 IR-mediated gene rearrangement (SCI Figure)
# Hypsizygus marmoreus mitochondrial gene arrangement diversity
# =============================================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(ggsci)
  library(cowplot)
})

# ── 1. 读取数据 ──────────────────────────────────────────────────────────────
genes <- read.csv("gene_orders.csv", stringsAsFactors = FALSE)

# ── 2. 颜色方案 (ggsci 风格) ─────────────────────────────────────────────────
gene_colors <- c(
  "cox1"  = "#0072B5FF",  # Lancet 蓝
  "cox2"  = "#20854EFF",  # Lancet 绿
  "cox3"  = "#E18727FF",  # Lancet 橙
  "nad1"  = "#BC3C29FF",  # Lancet 红
  "nad2"  = "#7876B1FF",  # Lancet 紫
  "nad3"  = "#6F99ADFF",  # Lancet 灰蓝
  "nad4"  = "#EE4C97FF",  # NPG 粉
  "nad4L" = "#FFDC91FF",  # 浅金
  "nad5"  = "#B09C85FF",  # JCO 棕
  "nad6"  = "#00A087FF",  # NPG 青绿
  "atp6"  = "#3C5488FF",  # NPG 藏蓝
  "atp8"  = "#F39B7FFF",  # NPG 珊瑚
  "atp9"  = "#91D1C2FF",  # NPG 浅青
  "cob"   = "#8491B4FF",  # NPG 灰紫
  "rps3"  = "#DC0000FF",  # NPG 红
  "dpo"   = "#B24745FF",  # 暗红
  "rnl"   = "#FF6F00",    # 亮橙 (突出)
  "rns"   = "#CE93D8",    # 浅紫
  "tRNA"  = "#B0C4DE"     # 浅钢蓝
)

# ── 3. 数据预处理 ────────────────────────────────────────────────────────────
n_patterns <- max(genes$pattern_id)
genes$y_pos <- n_patterns + 1 - genes$pattern_id
genes$x_pos <- genes$gene_idx
genes$color_key <- ifelse(genes$gene_type == "tRNA", "tRNA", genes$gene_name)

# ── 4. 构建反转区域 (动态检测 + 智能分类) ────────────────────────────────────

# 提取参考 (Pattern 1) 的基因链方向
ref <- genes %>% filter(pattern_id == 1, gene_type != "tRNA")
ref_strands <- setNames(ref$strand, ref$gene_name)

# 检测每个pattern中链方向翻转的CDS/rRNA基因(排除rnl单独处理)
flipped <- genes %>%
  filter(gene_type != "tRNA", gene_name != "rnl") %>%
  mutate(ref_strand = ref_strands[gene_name]) %>%
  filter(!is.na(ref_strand), strand != ref_strand)

# 将翻转基因分组: 找连续block (允许因tRNA间隔产生的gap <= 12)
find_blocks <- function(indices) {
  if (length(indices) == 0) return(data.frame(idx = integer(0), block = integer(0)))
  indices <- sort(indices)
  block_id <- cumsum(c(1, diff(indices) > 12))
  data.frame(idx = indices, block = block_id)
}

event_a_marker <- c("nad2", "nad3", "rps3")
event_b_marker <- c("atp8", "nad5", "nad4L", "cox3")

inv_regions <- data.frame()

for (pid in unique(flipped$pattern_id)) {
  sub <- flipped %>% filter(pattern_id == pid) %>% arrange(gene_idx)
  if (nrow(sub) == 0) next
  
  blocks <- find_blocks(sub$gene_idx)
  
  for (b in unique(blocks$block)) {
    block_idx <- blocks$idx[blocks$block == b]
    block_genes <- sub$gene_name[sub$gene_idx %in% block_idx]
    
    has_a <- any(block_genes %in% event_a_marker)
    has_b <- any(block_genes %in% event_b_marker)
    y <- n_patterns + 1 - pid
    
    if (has_a && has_b) {
      # 混合block: 拆分为 Event A 和 Event B 子区域
      a_idx <- sub$gene_idx[sub$gene_name %in% c(event_a_marker, "nad6", "atp6") & 
                            sub$gene_idx %in% block_idx]
      b_idx <- sub$gene_idx[sub$gene_name %in% c(event_b_marker, "dpo") & 
                            sub$gene_idx %in% block_idx]
      
      if (length(a_idx) > 0) {
        inv_regions <- rbind(inv_regions, data.frame(
          xmin = min(a_idx) - 0.08, xmax = max(a_idx) + 0.96,
          ymin = y - 0.46, ymax = y + 0.46,
          event = "Event A", pattern_id = pid
        ))
      }
      if (length(b_idx) > 0) {
        inv_regions <- rbind(inv_regions, data.frame(
          xmin = min(b_idx) - 0.08, xmax = max(b_idx) + 0.96,
          ymin = y - 0.46, ymax = y + 0.46,
          event = "Event B", pattern_id = pid
        ))
      }
    } else if (has_a) {
      inv_regions <- rbind(inv_regions, data.frame(
        xmin = min(block_idx) - 0.08, xmax = max(block_idx) + 0.96,
        ymin = y - 0.46, ymax = y + 0.46,
        event = "Event A", pattern_id = pid
      ))
    } else if (has_b) {
      inv_regions <- rbind(inv_regions, data.frame(
        xmin = min(block_idx) - 0.08, xmax = max(block_idx) + 0.96,
        ymin = y - 0.46, ymax = y + 0.46,
        event = "Event B", pattern_id = pid
      ))
    }
  }
}

# 扩展反转区域: 确保包含反转CDS基因之间的所有tRNA
for (i in seq_len(nrow(inv_regions))) {
  pid <- inv_regions$pattern_id[i]
  cds_min <- round(inv_regions$xmin[i] + 0.08)
  cds_max <- round(inv_regions$xmax[i] - 0.96)
  
  all_in_range <- genes %>% 
    filter(pattern_id == pid, gene_idx >= cds_min, gene_idx <= cds_max)
  
  if (nrow(all_in_range) > 0) {
    inv_regions$xmin[i] <- min(all_in_range$gene_idx) - 0.08
    inv_regions$xmax[i] <- max(all_in_range$gene_idx) + 0.96
  }
}

# rnl翻转标记
rnl_flip <- genes %>% 
  filter(gene_name == "rnl") %>%
  mutate(ref_strand = ref_strands["rnl"]) %>%
  filter(strand != ref_strand) %>%
  mutate(
    y_pos = n_patterns + 1 - pattern_id,
    xmin = gene_idx - 0.12,
    xmax = gene_idx + 1.04,
    ymin = y_pos - 0.50,
    ymax = y_pos + 0.50
  )

cat("检测到的反转区域:\n")
print(inv_regions %>% select(pattern_id, event, xmin, xmax))
cat("\nrnl翻转 (pattern, idx, strand):\n")
print(rnl_flip %>% select(pattern_id, gene_idx, strand))

# ── 5. 箭头多边形 ────────────────────────────────────────────────────────────

make_arrow_df <- function(x, y, strand, gene_type, w = 0.92, group_id) {
  h <- ifelse(gene_type == "tRNA", 0.25, 0.35)
  a <- 0.12
  
  if (strand == "+") {
    df <- data.frame(
      px = c(x, x + w - a, x + w, x + w - a, x),
      py = c(y + h, y + h, y, y - h, y - h)
    )
  } else {
    df <- data.frame(
      px = c(x + a, x + w, x + w, x + a, x),
      py = c(y + h, y + h, y - h, y - h, y)
    )
  }
  df$group_id <- group_id
  return(df)
}

arrow_list <- lapply(seq_len(nrow(genes)), function(i) {
  gid <- paste0("P", genes$pattern_id[i], "_G", genes$gene_idx[i])
  arr <- make_arrow_df(
    genes$x_pos[i], genes$y_pos[i], genes$strand[i],
    genes$gene_type[i], group_id = gid
  )
  arr$color_key <- genes$color_key[i]
  arr$gene_type <- genes$gene_type[i]
  arr
})
arrow_df <- do.call(rbind, arrow_list)

# ── 6. 标签数据 ───────────────────────────────────────────────────────────────

label_df <- genes %>%
  mutate(
    x_center = x_pos + 0.46,
    y_center = y_pos,
    display_label = gene_label
  )

# ── 7. Y轴标签 ──────────────────────────────────────────────────────────────

pattern_info <- data.frame(
  y = n_patterns:1,
  pid = 1:n_patterns,
  label = c(
    "Pattern 1: Reference (n=14)",
    "Pattern 2: Event A (n=11)",
    "Pattern 3: Event B + rnl(+) (n=2)",
    "Pattern 4: Event B variant (n=2)",
    "Pattern 5: Event A + B (n=1)",
    "Pattern 6: rnl flip only (n=1)",
    "Pattern 7: Ref \u2212 trnL/E (n=1)",
    "Pattern 8: Event A \u2212 trnL/E (n=1)"
  )
)

# ── 8. 图例顺序 ─────────────────────────────────────────────────────────────

legend_order <- c("cox1", "cox2", "cox3", "nad1", "nad2", "nad3", "nad4",
                  "nad4L", "nad5", "nad6", "atp6", "atp8", "atp9",
                  "cob", "rps3", "dpo", "rnl", "rns", "tRNA")

# ── 9. 主绑图 ────────────────────────────────────────────────────────────────

p <- ggplot() +
  
  # Event A 反转背景 (橙色虚线框)
  geom_rect(data = inv_regions %>% filter(event == "Event A"),
            aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
            fill = "#FFF8E1", color = "#E65100", linewidth = 0.55,
            linetype = "dashed", alpha = 0.7) +
  
  # Event B 反转背景 (蓝色虚线框)
  geom_rect(data = inv_regions %>% filter(event == "Event B"),
            aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
            fill = "#E8EAF6", color = "#1565C0", linewidth = 0.55,
            linetype = "dashed", alpha = 0.7) +
  
  # rnl 翻转高亮 (红色实线框)
  {if (nrow(rnl_flip) > 0)
    geom_rect(data = rnl_flip,
              aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
              fill = NA, color = "#C62828", linewidth = 1.0, linetype = "solid")
  } +
  
  # 基因箭头
  geom_polygon(data = arrow_df,
               aes(x = px, y = py, group = group_id, fill = color_key),
               color = "gray20", linewidth = 0.15) +
  
  # CDS/rRNA 基因标签
  geom_text(data = label_df %>% filter(gene_type != "tRNA"),
            aes(x = x_center, y = y_center, label = display_label),
            size = 1.9, color = "black") +
  
  # tRNA 标签
  geom_text(data = label_df %>% filter(gene_type == "tRNA"),
            aes(x = x_center, y = y_center, label = display_label),
            size = 1.4, color = "gray15") +
  
  # 颜色映射
  scale_fill_manual(
    values = gene_colors,
    breaks = legend_order,
    name = "Gene"
  ) +
  
  # 坐标轴
  scale_x_continuous(expand = expansion(mult = c(0.005, 0.005)), breaks = NULL) +
  scale_y_continuous(
    breaks = pattern_info$y,
    labels = pattern_info$label,
    expand = expansion(mult = c(0.05, 0.06))
  ) +
  
  # 标题
  labs(
    x = "Gene position (sorted by genomic location)",
    y = NULL,
    title = expression(paste(
      "Mitochondrial gene arrangement patterns in ",
      italic("Hypsizygus marmoreus")
    )),
    subtitle = "IR-mediated inversions generate 8 distinct gene orders across 31 genomes"
  ) +
  
  # 主题
  theme_minimal(base_size = 10) +
  theme(
    plot.title = element_text(size = 12, face = "bold", hjust = 0, margin = margin(b = 2)),
    plot.subtitle = element_text(size = 9, color = "gray30", hjust = 0, margin = margin(b = 8)),
    axis.text.y = element_text(size = 7, hjust = 1, color = "gray10"),
    axis.title.x = element_text(size = 8.5, margin = margin(t = 6)),
    axis.ticks = element_blank(),
    panel.grid = element_blank(),
    legend.position = "right",
    legend.title = element_text(size = 8.5, face = "bold"),
    legend.text = element_text(size = 7),
    legend.key.size = unit(0.38, "cm"),
    legend.spacing.y = unit(0.02, "cm"),
    plot.margin = margin(t = 8, r = 5, b = 8, l = 5),
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA)
  ) +
  
  guides(fill = guide_legend(ncol = 1,
                             override.aes = list(color = "gray30", linewidth = 0.2)))

# ── 10. 事件标注 ─────────────────────────────────────────────────────────────

# Event A 标注 (Pattern 2)
p2_a <- inv_regions %>% filter(pattern_id == 2, event == "Event A")
if (nrow(p2_a) > 0) {
  p <- p + annotate("text",
    x = (p2_a$xmin[1] + p2_a$xmax[1]) / 2,
    y = p2_a$ymax[1] + 0.10,
    label = "Event A  (18.5 kb, 385-bp IR)",
    size = 1.9, color = "#E65100", fontface = "bold")
}

# Event B 标注 (Pattern 3)
p3_b <- inv_regions %>% filter(pattern_id == 3, event == "Event B")
if (nrow(p3_b) > 0) {
  p <- p + annotate("text",
    x = (p3_b$xmin[1] + p3_b$xmax[1]) / 2,
    y = p3_b$ymax[1] + 0.10,
    label = "Event B  (12.2 kb, DR\u2192IR)",
    size = 1.9, color = "#1565C0", fontface = "bold")
}

# rnl flip 标注
if (nrow(rnl_flip) > 0) {
  for (i in seq_len(nrow(rnl_flip))) {
    p <- p + annotate("text",
      x = (rnl_flip$xmin[i] + rnl_flip$xmax[i]) / 2,
      y = rnl_flip$ymax[i] + 0.10,
      label = "rnl flip",
      size = 1.6, color = "#C62828", fontface = "bold.italic")
  }
}

# ── 11. 底部图注 ─────────────────────────────────────────────────────────────

p <- p + labs(caption = paste0(
  "Orange dashed box: Event A inversion (nad2-nad3-rps3 region, mediated by 385-bp inverted repeat, 88.1% identity).  ",
  "Blue dashed box: Event B inversion (atp8-nad5-nad4L-cox3 region, mediated by 213-bp DR-derived IR, 92.0%).\n",
  "Red solid box: rnl transcription direction flip from (\u2212) to (+).  ",
  "Arrow direction: \u2192 plus strand, \u2190 minus strand.  ",
  "Small blocks: tRNA genes (single-letter amino acid code)."
)) +
  theme(plot.caption = element_text(size = 6.5, color = "gray35", hjust = 0,
                                    lineheight = 1.4, margin = margin(t = 10)))

# ── 12. 保存 ─────────────────────────────────────────────────────────────────

ggsave("gene_rearrangement_figure.pdf", p,
       width = 36, height = 18, units = "cm", dpi = 600)
ggsave("gene_rearrangement_figure.png", p,
       width = 36, height = 18, units = "cm", dpi = 300)

cat("\n\u2713 Figure saved:\n")
cat("  gene_rearrangement_figure.pdf (36 \u00d7 18 cm, 600 dpi)\n")
cat("  gene_rearrangement_figure.png (36 \u00d7 18 cm, 300 dpi)\n")
