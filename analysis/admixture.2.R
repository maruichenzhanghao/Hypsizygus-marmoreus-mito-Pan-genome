# 加载必要的R包
library(tidyverse)   # 数据处理和绘图核心包
library(viridis)     # 科研级颜色方案
library(ggh4x)       # 精细分面控制（用于标签显示）

# 1. 设置工作目录（请替换为你的实际路径）
setwd("/home/maruichen/workspace/work/06.mito.sv/01.snp/02.sanjiake_wushaixuan/beifen1/")
# setwd("/home/maruichen/workspace/work/06.mito.sv/01.snp/01.sjianke/04.zhongqunjiegou1")


# 2. 读取样本ID（来自all.fam，保持原始顺序）
sample_info <- read_delim("all.fam", 
                          delim = " ", 
                          col_names = c("FamID", "IndID", "Paternal", "Maternal", "Sex", "Pheno"),
                          show_col_types = FALSE)
sample_ids <- sample_info$IndID  # 提取样本ID，用于后续标注


# 3. 批量读取K=2到K=8的Q文件并整理数据
k_values <- 2:8                  # 定义K值范围
max_k <- max(k_values)           # 最大K值（最后一层显示标签）
admixture_data <- list()         # 存储所有K值的数据

for (k in k_values) {
  # 读取Q文件（空格分隔，无表头）
  q_file <- paste0("all.", k, ".Q")
  q_data <- read.table(q_file, header = FALSE, sep = " ")
  
  # 为列命名（Pop1到Popk）
  colnames(q_data) <- paste0("Pop", 1:k)
  
  # 添加样本ID和K值信息
  q_data <- q_data %>%
    mutate(
      Sample = sample_ids,  # 匹配样本ID（保持原始顺序）
      K = k,                # 标记当前K值
      Sample = factor(Sample, levels = sample_ids)  # 锁定样本顺序
    ) %>%
    # 转换为长格式（适合ggplot绘图）
    pivot_longer(
      cols = starts_with("Pop"),  # 选择所有群体列
      names_to = "Population",    # 群体名称（Pop1-Popk）
      values_to = "Proportion"    # 祖先成分比例
    )
  
  # 存入列表
  admixture_data[[as.character(k)]] <- q_data
}

# 合并所有K值的数据
all_data <- bind_rows(admixture_data)


# 4. 定义固定颜色映射（确保同一群体颜色统一）
unique_pops <- sort(unique(all_data$Population))  # 提取所有群体并排序（Pop1到Pop8）
pop_colors <- viridis(length(unique_pops), option = "turbo")  # 分配颜色（数量与群体一致）


# 5. 绘制Admixture图（新增边框设置）
admixture_plot <- ggplot(all_data, aes(x = Sample, y = Proportion, fill = Population)) +
  # 堆叠条形图 + 边框设置 ↓↓↓
  geom_col(
    position = "fill", 
    width = 1, 
    color = "black",    # 边框颜色（可改为#333333更柔和）
    size = 0.2          # 线条粗细（0.2~0.5之间调整）
  ) +
  
  # 按K值分面（纵向排列）
  facet_wrap(~K, ncol = 1, scales = "free_x") +
  
  # 固定颜色映射（同一群体颜色统一）
  scale_fill_manual(
    values = pop_colors,
    breaks = unique_pops,
    name = "Ancestry"  # 图例标题
  ) +
  
  # 仅在最后一层（最大K值）显示样本标签
  facetted_pos_scales(
    x = list(
      K == max_k ~ scale_x_discrete(labels = sample_ids),  # 最后一层显示标签
      K != max_k ~ scale_x_discrete(labels = NULL)         # 其他层隐藏标签
    )
  ) +
  
  # 坐标轴和标题设置
  labs(
    title = "Admixture Ancestry Components Across Different K Values",
    x = "Sample",
    y = "Proportion"
  ) +
  
  # 主题美化
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold", margin = margin(b = 15)),
    axis.title.x = element_text(size = 12, face = "bold", margin = margin(t = 10)),
    axis.title.y = element_text(size = 12, face = "bold", margin = margin(r = 10)),
    # 样本标签旋转（仅最后一层有效）
    axis.text.x = element_text(
      angle = 45, 
      hjust = 1, 
      vjust = 0.5, 
      size = 7  # 样本多，缩小字体避免重叠
    ),
    strip.text.x = element_text(size = 11, face = "bold", color = "#333333"),
    legend.position = "right",
    legend.title = element_text(size = 11, face = "bold"),
    legend.text = element_text(size = 9),
    plot.margin = margin(10, 10, 20, 10)
  )

# 显示图形
print(admixture_plot)
 
# 保存图形（PDF矢量格式）
ggsave(
  "admixture_all_K.shaixuan2.pdf",
  plot = admixture_plot,
  width = 12,  # 宽度根据样本数量调整
  height = 2.5 * length(k_values),  # 每个K值分配2.5英寸高度
  dpi = 600,
  device = "pdf"
)

ggsave(
  "admixture_all_K.shaixuan1.5.pdf",
  plot = admixture_plot,
  width = 12,  # 宽度根据样本数量调整
  height = 1.5 * length(k_values),  # 每个K值分配2.5英寸高度
  dpi = 600,
  device = "pdf"
)

