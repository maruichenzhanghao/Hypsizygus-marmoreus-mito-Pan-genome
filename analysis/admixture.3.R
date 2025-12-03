# 加载必要的R包
library(tidyverse)   # 数据处理和绘图核心包
library(viridis)     # 科研级颜色方案
library(ggh4x)       # 精细分面控制

# 1. 设置工作目录（请替换为你的实际路径）
setwd("/home/maruichen/workspace/work/06.mito.sv/01.snp/02.sanjiake_wushaixuan/beifen1/")


# 2. 读取群体信息文件（pop.txt）
pop_info <- read_delim(
  "../pop.txt",
  delim = " ",  # 按空格分隔
  col_names = c("Sample", "Group"),  # 定义列名：样本ID和群体编号
  show_col_types = FALSE
)

# 按群体编号（1-5）排序，得到最终样本顺序
sorted_samples <- pop_info %>%
  arrange(Group) %>%  # 按Group列（1-5）排序
  pull(Sample)  # 提取排序后的样本ID列表


# 3. 读取样本原始ID（用于匹配Q文件）
sample_info <- read_delim(
  "all.fam", 
  delim = " ", 
  col_names = c("FamID", "IndID", "Paternal", "Maternal", "Sex", "Pheno"),
  show_col_types = FALSE
)
sample_ids <- sample_info$IndID  # 原始样本ID


# 4. 批量读取K=2到K=8的Q文件并按群体排序
k_values <- 2:8
admixture_data <- map_dfr(k_values, function(k) {
  q_file <- paste0("all.", k, ".Q")
  q_data <- read.table(q_file, header = FALSE, sep = " ")
  colnames(q_data) <- paste0("Pop", 1:k)  # 命名群体列
  
  q_data %>% 
    mutate(Sample = sample_ids) %>%  # 关联样本ID
    # 按pop.txt中的群体顺序排序
    arrange(factor(Sample, levels = sorted_samples)) %>%
    pivot_longer(
      cols = -Sample,
      names_to = "Population",
      values_to = "Proportion"
    ) %>%
    mutate(
      Sample = factor(Sample, levels = sorted_samples),  # 锁定排序
      K = k
    )
})


# 5. 定义自定义颜色方案
# 获取所有群体并排序（确保顺序一致）
unique_pops <- sort(unique(admixture_data$Population))  

# 步骤1：预定义1-5号群体的颜色（命名向量）
custom_colors <- c(
  "Pop1" = "#BCE266",  # 1号群体颜色
  "Pop2" = "#FF6767",  # 2号群体颜色
  "Pop3" = "#FFBF66",  # 3号群体颜色
  "Pop4" = "#ACCCEB",  # 4号群体颜色
  "Pop5" = "#AA9CCD"   # 5号群体颜色
)

# 步骤2：筛选需要用viridis配色的群体（6号及以上）
viridis_pops <- unique_pops[!unique_pops %in% names(custom_colors)]

# 步骤3：为这些群体生成viridis颜色（数量严格匹配）
viridis_colors <- viridis(length(viridis_pops), option = "turbo")  
names(viridis_colors) <- viridis_pops  # 命名，确保一一对应

# 步骤4：合并颜色，保持unique_pops的顺序
pop_colors <- sapply(unique_pops, function(pop) {
  if (pop %in% names(custom_colors)) {
    custom_colors[pop]  # 优先用自定义颜色
  } else {
    viridis_colors[pop]  # 其余用viridis
  }
})
# 确保颜色向量顺序与群体顺序一致（可选，但若群体命名规则一致可省略）
names(pop_colors) <- unique_pops


# 6. 绘制Admixture图（无分隔线+自定义颜色）
admixture_plot <- ggplot(admixture_data, aes(x = Sample, y = Proportion, fill = Population)) +
  geom_col(
    position = "fill", 
    width = 1, 
    color = "black", 
    size = 0.2
  ) +
  
  # 按K值分面
  facet_wrap(~K, ncol = 1, scales = "free_x") +
  
  # 应用自定义颜色
  scale_fill_manual(
    values = pop_colors,
    breaks = unique_pops,
    name = "Ancestry"
  ) +
  
  # 仅在最大K值显示样本标签
  facetted_pos_scales(
    x = list(
      K == max(k_values) ~ scale_x_discrete(labels = sorted_samples),
      TRUE ~ scale_x_discrete(labels = NULL)
    )
  ) +
  
  labs(
    title = "Admixture Ancestry (Grouped by pop.txt)",
    x = "Sample",
    y = "Proportion"
  ) +
  
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 7),
    strip.text.x = element_text(size = 11, face = "bold"),
    legend.position = "right"
  )

print(admixture_plot)

# 保存图形（可选）
ggsave("admixture_custom_colors.png", width = 14, height = 10, dpi = 300)

ggsave(
  "admixture_all1.5.pdf",
  plot = admixture_plot,
  width = 12,  # 宽度根据样本数量调整
  height = 1.5 * length(k_values),  # 每个K值分配2.5英寸高度
  dpi = 600,
  device = "pdf"
)
