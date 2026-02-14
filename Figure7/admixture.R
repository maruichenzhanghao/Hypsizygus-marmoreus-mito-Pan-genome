# 加载必要的R包
library(tidyverse)   # 包含dplyr、ggplot2、tidyr等，用于数据处理和绘图
library(viridis)     # 科研级颜色方案
library(patchwork)   # 用于组合多个图形（可选）

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
admixture_data <- list()  # 存储所有K值的数据

for (k in 2:8) {
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


# 4. 绘制Admixture图（按K值分面）
admixture_plot <- ggplot(all_data, aes(x = Sample, y = Proportion, fill = Population)) +
  # 堆叠条形图（position="fill"确保每个样本的总比例为1）
  geom_col(position = "fill", width = 1) +
  
  # 按K值分面（每个K值一个子图）
  facet_wrap(~K, ncol = 1, scales = "free_x") +  # ncol=1表示纵向排列，可改为nrow=2横向
  
  # 颜色方案（使用viridis，区分度高且色盲友好）
  scale_fill_viridis(discrete = TRUE, option = "turbo", name = "Ancestry") +
  
  # 坐标轴和标题设置
  labs(
    title = "Admixture Ancestry Components Across Different K Values",
    x = "Sample",
    y = "Proportion"
  ) +
  
  # 主题美化
  theme_minimal() +
  theme(
    # 标题居中
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold", margin = margin(b = 15)),
    # 坐标轴标签
    axis.title.x = element_text(size = 12, face = "bold", margin = margin(t = 10)),
    axis.title.y = element_text(size = 12, face = "bold", margin = margin(r = 10)),
    # 样本标签旋转（避免重叠）
    axis.text.x = element_text(
      angle = 90, 
      hjust = 1, 
      vjust = 0.5, 
      size = 7  # 样本多，缩小字体
    ),
    # 分面标题
    strip.text.x = element_text(size = 11, face = "bold", color = "#333333"),
    # 图例
    legend.position = "right",
    legend.title = element_text(size = 11, face = "bold"),
    legend.text = element_text(size = 9),
    # 调整边距
    plot.margin = margin(10, 10, 20, 10)
  )

# 显示图形
print(admixture_plot)

# 保存图形（PDF格式支持矢量图，方便后续编辑）
ggsave(
  "admixture_all_K.shaixuan.pdf",
  plot = admixture_plot,
  width = 12,  # 宽度根据样本数量调整（样本多则加宽）
  height = 2.5 * length(2:8),  # 每个K值分配2.5英寸高度
  dpi = 300,
  device = "pdf"
)
