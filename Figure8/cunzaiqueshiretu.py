import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import glob
import os

# 找到最新的OrthoFinder结果目录
results_dir = max(glob.glob(os.path.expanduser("~/workspace/01.mito/02.gb/fasta_output/OrthoFinder/Results_*")), key=os.path.getctime)
gene_count_file = f"{results_dir}/Orthogroups/Orthogroups.GeneCount.tsv"

# 读取基因计数文件
df = pd.read_csv(gene_count_file, sep="\t")

# 创建二进制矩阵（1=存在，0=缺失）
binary_matrix = (df.iloc[:, 1:-1] > 0).astype(int)
binary_matrix.index = df["Orthogroup"]

# 过滤掉全0或全1的orthogroups（聚焦壳基因）
binary_matrix = binary_matrix[(binary_matrix.sum(axis=1) > 0) & (binary_matrix.sum(axis=1) < 31)]

# 绘制热图
plt.figure(figsize=(12, 8))
sns.heatmap(binary_matrix, cmap="Blues", cbar=False)
plt.xlabel("Strains")
plt.ylabel("Orthogroups")
plt.title("Presence/Absence of Shell Genes in Lentinula edodes mtDNA")
plt.tight_layout()
plt.savefig(os.path.expanduser("~/workspace/01.mito/02.gb/pangenome_heatmap.png"))
plt.show()
