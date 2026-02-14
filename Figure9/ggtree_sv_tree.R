#!/usr/bin/env Rscript
# ============================================================
# 真姬菇线粒体基因组 SV-based 系统发育树 (ggtree 科研绘图)
# Based on Structural Variation genotype data (62 SV sites)
# ============================================================

library(ggtree)
library(ggplot2)
library(treeio)
library(ape)

# --- 1. 读取Newick树 ---
tree <- read.tree("sv_nj_tree.nwk")

# --- 2. 样本地理来源元数据 ---
meta <- data.frame(
  label = c(
    # 用户自有样品 - 上海 (f2是参考不在树中)
    "f4", "nn12-1", "nn12-17",
    # GenBank参考基因组 - 福州
    "MF133443", "MH382825",
    # NCBI SRA - 韩国
    "SRR7874787",
    # NCBI SRA - 福建
    "SRR12151860", "SRR12151871", "SRR12151875", "SRR12151883",
    # NCBI SRA - 福州 (FQX_MS01) - 18个样品
    "SRR8699796", "SRR8699797", "SRR8699800", "SRR8699801",
    "SRR8699802", "SRR8699803", "SRR8699804", "SRR8699805",
    "SRR8699808", "SRR8699809", "SRR8699811", "SRR8699813",
    "SRR8699814", "SRR8699815", "SRR8699816", "SRR8699817",
    "SRR8699833", "SRR8699834", "SRR8699835", "SRR8699837"
  ),
  stringsAsFactors = FALSE
)

# 地理来源
meta$Origin <- ifelse(meta$label %in% c("f4", "nn12-1", "nn12-17"), "Shanghai",
               ifelse(meta$label %in% c("MF133443", "MH382825"), "Fuzhou",
               ifelse(meta$label == "SRR7874787", "Korea",
               ifelse(meta$label %in% c("SRR12151860","SRR12151871","SRR12151875","SRR12151883"), "Fujian",
                      "Fuzhou"))))

# 数据来源
meta$Source <- ifelse(meta$label %in% c("f4", "nn12-1", "nn12-17"), "This study",
               ifelse(meta$label %in% c("MF133443", "MH382825"), "GenBank",
                      "NCBI SRA"))

# 合并福建/福州为一个大区域类别用于颜色标注
meta$Region <- ifelse(meta$Origin == "Shanghai", "Shanghai, China",
               ifelse(meta$Origin == "Korea", "Seoul, Korea",
                      "Fujian, China"))

# --- 3. f2是参考，不在树里，其余30个样本 ---
# 注意: f2在构建树时是参考路径，不包含在VCF样本中
# f4 在树中

# --- 4. 颜色方案 (科研配色) ---
region_colors <- c(
  "Shanghai, China" = "#E64B35",   # 红色 - Nature配色
  "Fujian, China"   = "#4DBBD5",   # 蓝色
  "Seoul, Korea"    = "#00A087"    # 绿色
)

source_shapes <- c(
  "This study" = 16,   # 实心圆
  "GenBank"    = 17,   # 三角形
  "NCBI SRA"   = 15    # 方形
)

# --- 5. 绑定元数据到树 ---
p_base <- ggtree(tree, layout = "rectangular", 
                 size = 0.8, color = "black") %<+% meta

# --- 6. 绘图 ---
p <- p_base +
  # 树枝末端标注 (彩色点，按地区上色，按来源设形状)
  geom_tippoint(aes(color = Region, shape = Source), size = 3.5, stroke = 0.5) +
  
  # 样本标签
  geom_tiplab(aes(color = Region), 
              size = 3.2, hjust = -0.15, 
              fontface = "plain") +
  
  # 标尺
  geom_treescale(x = 0, y = -1.5, width = 0.05, 
                 offset = 0.5, fontsize = 3.5,
                 linesize = 0.8) +
  
  # 颜色和形状
  scale_color_manual(values = region_colors, name = "Geographic origin") +
  scale_shape_manual(values = source_shapes, name = "Data source") +
  
  # 主题
  theme_tree2() +
  theme(
    # 图例
    legend.position = c(0.25, 0.85),
    legend.background = element_rect(fill = "white", color = "grey80", linewidth = 0.3),
    legend.title = element_text(size = 11, face = "bold"),
    legend.text = element_text(size = 10),
    legend.key.size = unit(0.5, "cm"),
    legend.spacing.y = unit(0.1, "cm"),
    # 坐标轴
    axis.text.x = element_text(size = 10),
    # 整体
    plot.margin = margin(10, 80, 10, 10),
    plot.title = element_text(size = 13, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 10, hjust = 0.5, color = "grey40")
  ) +
  
  # 标题
  labs(
    title = expression(italic("Hypsizygus marmoreus") ~ "mitochondrial genome phylogeny"),
    subtitle = "Neighbor-joining tree based on 62 structural variation sites"
  ) +
  
  # 留出标签空间
  xlim(NA, max(fortify(tree)$x) * 1.8)

# --- 7. 输出PDF和PNG ---
ggsave("sv_nj_tree_ggtree.pdf", p, width = 10, height = 10, dpi = 300)
ggsave("sv_nj_tree_ggtree.png", p, width = 10, height = 10, dpi = 300)

cat("Done! Output files:\n")
cat("  sv_nj_tree_ggtree.pdf\n")
cat("  sv_nj_tree_ggtree.png\n")
