library(plotly)
# install.packages("plotly")
setwd("/home/maruichen/workspace/work/06.mito.sv/01.snp/02.sanjiake_wushaixuan/")

# 读取 PCA 结果
pca <- read.table("outputPrefix.eigenvec", header = FALSE)
colnames(pca) <- c("FID", "IID", paste0("PC", 1:(ncol(pca) - 2)))

# 读取群体信息
pop <- read.table("pop.txt", header = FALSE)
colnames(pop) <- c("IID", "Population")

# 合并
df <- merge(pca, pop, by = "IID")

# 画 3D 图
fig <- plot_ly(
  df,
  x = ~PC1, y = ~PC2, z = ~PC3,
  color = ~as.factor(Population),
  colors = "Set1",
  text = ~IID,
  type = "scatter3d",
  mode = "markers"
) %>%
  layout(scene = list(
    xaxis = list(title = "PC1"),
    yaxis = list(title = "PC2"),
    zaxis = list(title = "PC3")
  ))

fig



library(rgl)
# install.packages("rgl")


# 绘制
plot3d(df$PC1, df$PC2, df$PC3,
       col = as.numeric(as.factor(df$Population)),
       size = 5, type = "s")
legend3d("topright", legend = unique(df$Population),
         col = 1:length(unique(df$Population)), pch = 16)
