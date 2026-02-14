import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import glob
import os

# 找到最新的OrthoFinder结果目录
results_dir = max(glob.glob(os.path.expanduser("~/workspace/01.mito/02.gb/fasta_output/OrthoFinder/Results_*")), key=os.path.getctime)
gene_count_file = f"{results_dir}/Orthogroups/Orthogroups.GeneCount.tsv"

# 读取基因计数文件
df = pd.read_csv(gene_count_file, sep="\t")

# 计算每个株系的核心和壳基因数量
strain_counts = (df.iloc[:, 1:-1] > 0).sum(axis=0)  # 每个株系的总基因数
core_counts = (df.iloc[:, 1:-1] > 0).sum(axis=1)
core_mask = core_counts == 31  # 核心基因（存在于所有31株）
shell_mask = (core_counts >= 2) & (core_counts < 31)  # 壳基因

# 统计每个株系的核心和壳基因
core_per_strain = (df.loc[core_mask, df.columns[1:-1]] > 0).sum(axis=0)
shell_per_strain = (df.loc[shell_mask, df.columns[1:-1]] > 0).sum(axis=0)

# 创建数据框用于绘图
plot_data = pd.DataFrame({
    'Core': core_per_strain,
    'Shell': shell_per_strain
}, index=df.columns[1:-1])

# 绘制堆叠柱状图
plot_data.plot(kind='bar', stacked=True, color=['#1f77b4', '#ff7f0e'], figsize=(12, 6))
plt.xlabel("Strains")
plt.ylabel("Number of Orthogroups")
plt.title("Core and Shell Genes per Strain in Lentinula edodes mtDNA")
plt.xticks(rotation=45, ha='right')
plt.tight_layout()
plt.savefig(os.path.expanduser("~/workspace/01.mito/02.gb/stacked_pangenome.png"))
plt.show()
