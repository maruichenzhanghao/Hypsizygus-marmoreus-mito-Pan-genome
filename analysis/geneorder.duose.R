setwd("/home/maruichen/workspace/9.R/02.mito/08.geneorder")

# 加载所需包（确保优先加载dplyr，减少函数冲突）
library(ggplot2)
library(dplyr)       # 显式加载，优先处理数据操作
library(scales)      # 用于生成多色渐变
library(magrittr)    # 确保管道符%>%正常工作（部分环境依赖）

# 1. 读取预处理数据（添加数据读取检查，避免后续报错）
df <- read.csv("gene_order_for_R.csv", stringsAsFactors = FALSE)
if (!nrow(df) > 0) stop("数据读取失败！请检查文件路径或文件格式")  # 防错机制

# 2. 细分基因类别（区分tRNA、rRNA、蛋白编码基因）
# 显式使用dplyr::mutate和dplyr::case_when，彻底避免函数冲突
df <- df %>%
  dplyr::mutate(
    gene_category = dplyr::case_when(
      gene_type == "tRNA" ~ "tRNA",
      gene_name == "rnl" ~ "rRNA (rnl)",
      gene_name == "rns" ~ "rRNA (rns)",
      TRUE ~ "Protein-coding gene"  # 16个蛋白编码基因
    )
  )

# 3. 定义颜色方案（优化颜色映射逻辑，确保无遗漏）
# 提取蛋白编码基因列表，避免重复
protein_genes <- df %>%
  dplyr::filter(gene_category == "Protein-coding gene") %>%
  dplyr::pull(gene_name) %>%
  unique()

# 生成蛋白编码基因的唯一颜色（hue色板兼容性强，自动匹配数量）
protein_colors <- hue_pal()(length(protein_genes))
names(protein_colors) <- protein_genes  # 绑定基因名与颜色，避免错位

# 合并所有基因类别的颜色（按tRNA→rRNA→蛋白基因顺序，与图例一致）
color_palette <- c(
  "tRNA" = "#a8d1ff",          # tRNA：浅蓝色（保留原色）
  "rRNA (rnl)" = "#7ccc7c",    # rRNA (rnl)：浅绿色
  "rRNA (rns)" = "#66c2a5",    # rRNA (rns)：深绿色（与rnl区分更明显）
  protein_colors               # 蛋白编码基因：自动生成的多色
)

# 4. 计算紧凑x轴位置（按order_id分组，避免不同组位置重叠）
compact_factor <- 0.7  # 保留原紧凑系数，平衡显示密度
df <- df %>%
  dplyr::group_by(order_id) %>%
  dplyr::mutate(compact_pos = position * compact_factor) %>%
  dplyr::ungroup()  # 取消分组，避免后续操作受影响

# 5. 计算x轴边界（添加微小缓冲，避免基因块超出画布）
min_pos <- min(df$compact_pos) - 0.1
max_pos <- max(df$compact_pos) + 0.1

# 6. 提取基因组标注信息（显式指定dplyr::select，解决原冲突核心问题）
order_labels <- df %>%
  dplyr::select(order_id, genomes) %>%  # 明确用dplyr的select，不与其他包冲突
  unique() %>%
  dplyr::mutate(
    label_y = order_id * 2 - 1.2,       # 标注在每行下方，位置更合理
    label_x = (min_pos + max_pos) / 2   # 标注水平居中
  )

# 7. 创建绘图对象（优化图层顺序和图例逻辑）
p <- ggplot() +
  # 7.1 绘制蛋白编码基因（最底层，避免被覆盖）
  geom_rect(
    data = df %>% dplyr::filter(gene_category == "Protein-coding gene"),
    aes(
      xmin = compact_pos - 0.35, xmax = compact_pos + 0.35,
      ymin = order_id * 2 - 0.4, ymax = order_id * 2 + 0.4,
      fill = gene_name  # 按基因名着色（每个基因唯一颜色）
    ),
    color = "black", linewidth = 0.3  # 黑色边框，清晰区分基因块
  ) +
  # 7.2 绘制rRNA（中层，覆盖蛋白基因但不覆盖tRNA标注）
  geom_rect(
    data = df %>% dplyr::filter(gene_category %in% c("rRNA (rnl)", "rRNA (rns)")),
    aes(
      xmin = compact_pos - 0.35, xmax = compact_pos + 0.35,
      ymin = order_id * 2 - 0.4, ymax = order_id * 2 + 0.4,
      fill = gene_category  # 按rRNA类型着色
    ),
    color = "black", linewidth = 0.3
  ) +
  # 7.3 绘制tRNA（顶层，确保标注清晰）
  geom_rect(
    data = df %>% dplyr::filter(gene_category == "tRNA"),
    aes(
      xmin = compact_pos - 0.3, xmax = compact_pos + 0.3,  # tRNA宽度稍窄，区分类型
      ymin = order_id * 2 - 0.4, ymax = order_id * 2 + 0.4,
      fill = "tRNA"
    ),
    color = "black", linewidth = 0.3
  ) +
  # 7.4 添加基因标签（优化字体大小，避免重叠）
  # tRNA标签（显示trna_code，字体稍大）
  geom_text(
    data = df %>% dplyr::filter(gene_category == "tRNA"),
    aes(x = compact_pos, y = order_id * 2, label = trna_code),
    size = 2.2, color = "#333333"  # 深灰色标签，更醒目
  ) +
  # 非tRNA标签（显示gene_name，字体稍小）
  geom_text(
    data = df %>% dplyr::filter(gene_category != "tRNA"),
    aes(x = compact_pos, y = order_id * 2, label = gene_name),
    size = 1.8, color = "#333333"
  ) +
  # 7.5 添加基因组名称标注（避免与基因块重叠）
  geom_text(
    data = order_labels,
    aes(x = label_x, y = label_y, label = paste0("Genomes: ", genomes)),
    size = 2.5, hjust = 0.5, color = "#666666"  # 灰色标注，不抢焦点
  ) +
  # 7.6 坐标轴设置（优化y轴反转逻辑，确保order_id=1在最上方）
  scale_x_continuous(
    limits = c(min_pos, max_pos),
    expand = expansion(add = 0),  # 无额外扩展，避免空白
    name = "Gene Position (sorted by location)"  # x轴标题
  ) +
  scale_y_reverse(
    breaks = seq(2, max(df$order_id)*2, by = 2),  # y轴刻度与order_id对应
    labels = 1:max(df$order_id),                  # 显示1~n的序号
    expand = expansion(add = c(1, 2)),            # 上下留空，避免标注被切
    name = "Unique Gene Orders"                   # y轴标题
  ) +
  # 7.7 颜色映射（关键优化：确保图例与颜色方案完全匹配）
  scale_fill_manual(
    values = color_palette,
    name = "Gene Type",  # 图例标题
    breaks = c("tRNA", "rRNA (rnl)", "rRNA (rns)", protein_genes),  # 图例顺序
    drop = FALSE  # 强制显示所有颜色（避免部分基因不显示在图例）
  ) +
  # 7.8 主题设置（优化可读性和美观度）
  theme_minimal() +
  theme(
    # 隐藏x轴文字和刻度（x轴仅示意位置，无需具体数值）
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    # 隐藏网格线（避免干扰基因块显示）
    panel.grid = element_blank(),
    # 图例优化（右侧放置，大小适配）
    legend.position = "right",
    legend.key.size = unit(0.5, "cm"),    # 图例图标大小
    legend.text = element_text(size = 7), # 图例文字大小（避免拥挤）
    legend.title = element_text(size = 8),# 图例标题大小
    # 调整边距（右侧留足图例空间，底部留标注空间）
    plot.margin = margin(10, 80, 30, 10), # 上、右、下、左
    # 坐标轴标题大小优化
    axis.title.x = element_text(size = 9, margin = margin(t = 10)),
    axis.title.y = element_text(size = 9, margin = margin(r = 10))
  )

print(p)
# 8. 保存图形（按你说的“2比较合适”，保留宽12、高6的尺寸，优化分辨率）
# 只保留适配的版本（避免冗余文件，若需原版本可注释此句）
ggsave(
  filename = "gene_order_colored_final.svg",  # 重命名为“最终版”，便于识别
  plot = p,
  width = 12,    # 你认为合适的宽度
  height = 6,    # 你认为合适的高度
  device = "svg",# 矢量图格式，放大不失真
  dpi = 300      # 高分辨率，适配期刊/报告
)

# 9. 提示完成
cat("绘图完成！已生成优化版矢量图：gene_order_colored_final.svg\n")
