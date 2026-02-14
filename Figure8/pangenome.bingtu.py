import pandas as pd
import matplotlib.pyplot as plt
import glob
import os

# 找到最新的OrthoFinder结果目录
results_dir = max(glob.glob(os.path.expanduser("~/workspace/01.mito/02.gb/fasta_output/OrthoFinder/Results_*")), key=os.path.getctime)
gene_count_file = f"{results_dir}/Orthogroups/Orthogroups.GeneCount.tsv"

# 读取基因计数文件
df = pd.read_csv(gene_count_file, sep="\t")

# 计算核心和壳基因
core_counts = (df.iloc[:, 1:-1] > 0).sum(axis=1)
core = sum(core_counts == 31)
shell = sum((core_counts >= 2) & (core_counts < 31))

# 绘制饼图
labels = ["Core", "Shell"]
sizes = [core, shell]
colors = ["#1f77b4", "#ff7f0e"]
plt.pie(sizes, labels=labels, colors=colors, autopct='%1.1f%%', startangle=90)
plt.title("Core and Shell Genes in Lentinula edodes mtDNA")
plt.axis('equal')
plt.savefig(os.path.expanduser("~/workspace/01.mito/02.gb/pangenome_pie.png"))
plt.show()
