# 加载必要的包
library(ggplot2)
library(dplyr)

# 创建示例数据（线粒体密码子RSCU值）
codon_data <- data.frame(
  AminoAcid = rep(c("Ala", "Arg", "Asn", "Asp", "Cys", "Gln", "Glu", "Gly", "His", "Ile", "Leu", "Lys", "Met", "Phe", "Pro", "Ser", "Thr", "Trp", "Tyr", "Val", "Stop"), each = 6),
  Codon = c(
    "GCU", "GCC", "GCA", "GCG", NA, NA,  # Ala
    "CGU", "CGC", "CGA", "CGG", "AGA", "AGG",  # Arg
    "AAU", "AAC", NA, NA, NA, NA,  # Asn
    "GAU", "GAC", NA, NA, NA, NA,  # Asp
    "UGU", "UGC", NA, NA, NA, NA,  # Cys
    "CAA", "CAG", NA, NA, NA, NA,  # Gln
    "GAA", "GAG", NA, NA, NA, NA,  # Glu
    "GGU", "GGC", "GGA", "GGG", NA, NA,  # Gly
    "CAU", "CAC", NA, NA, NA, NA,  # His
    "AUU", "AUC", "AUA", NA, NA, NA,  # Ile
    "UUA", "UUG", "CUU", "CUC", "CUA", "CUG",  # Leu
    "AAA", "AAG", NA, NA, NA, NA,  # Lys
    "AUG", NA, NA, NA, NA, NA,  # Met
    "UUU", "UUC", NA, NA, NA, NA,  # Phe
    "CCU", "CCC", "CCA", "CCG", NA, NA,  # Pro
    "UCU", "UCC", "UCA", "UCG", "AGU", "AGC",  # Ser
    "ACU", "ACC", "ACA", "ACG", NA, NA,  # Thr
    "UGG", NA, NA, NA, NA, NA,  # Trp
    "UAU", "UAC", NA, NA, NA, NA,  # Tyr
    "GUU", "GUC", "GUA", "GUG", NA, NA,  # Val
    "UAA", "UAG", "UGA", NA, NA, NA  # Stop
  ),
  RSCU = c(
    0.8, 1.5, 1.2, 0.5, NA, NA,  # Ala
    0.7, 1.2, 0.9, 0.8, 1.5, 1.9,  # Arg
    0.9, 1.1, NA, NA, NA, NA,  # Asn
    0.8, 1.2, NA, NA, NA, NA,  # Asp
    0.7, 1.3, NA, NA, NA, NA,  # Cys
    0.6, 1.4, NA, NA, NA, NA,  # Gln
    0.7, 1.3, NA, NA, NA, NA,  # Glu
    0.9, 1.5, 1.2, 0.4, NA, NA,  # Gly
    0.8, 1.2, NA, NA, NA, NA,  # His
    1.1, 1.8, 0.1, NA, NA, NA,  # Ile
    0.4, 0.3, 0.7, 1.2, 0.8, 1.6,  # Leu
    0.7, 1.3, NA, NA, NA, NA,  # Lys
    1.0, NA, NA, NA, NA, NA,  # Met
    1.2, 0.8, NA, NA, NA, NA,  # Phe
    0.9, 1.3, 1.1, 0.7, NA, NA,  # Pro
    0.8, 1.3, 0.9, 0.5, 1.0, 1.5,  # Ser
    0.9, 1.4, 1.1, 0.6, NA, NA,  # Thr
    1.0, NA, NA, NA, NA, NA,  # Trp
    1.1, 0.9, NA, NA, NA, NA,  # Tyr
    0.9, 1.4, 1.0, 0.7, NA, NA,  # Val
    0.8, 0.2, 1.0, NA, NA, NA  # Stop
  )
)

# 移除NA值
codon_data <- na.omit(codon_data)

# 定义氨基酸顺序（按字母排序）
aa_order <- c("Ala", "Arg", "Asn", "Asp", "Cys", "Gln", "Glu", "Gly", 
              "His", "Ile", "Leu", "Lys", "Met", "Phe", "Pro", "Ser", 
              "Thr", "Trp", "Tyr", "Val", "Stop")

# 转换为因子以控制顺序
codon_data$AminoAcid <- factor(codon_data$AminoAcid, levels = aa_order)

# 绘制堆积柱状图
ggplot(codon_data, aes(x = AminoAcid, y = RSCU, fill = Codon)) +
  geom_bar(stat = "identity", color = "black", size = 0.2) +
  labs(
    title = "线粒体密码子RSCU值堆积柱状图",
    x = "氨基酸",
    y = "相对同义密码子使用度 (RSCU)",
    fill = "密码子"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
    legend.position = "bottom",
    legend.text = element_text(size = 8),
    legend.title = element_text(size = 10),
    panel.grid.major.x = element_blank(),
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold")
  ) +
  guides(fill = guide_legend(nrow = 7, byrow = TRUE))
