import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import glob
import os
from random import shuffle
import seaborn as sns  # 用于美化图形

# 设置Seaborn风格（可选：'whitegrid', 'darkgrid', 'ticks'）
sns.set_style("whitegrid")
plt.rcParams["font.family"] = "Arial"  # 推荐Arial或Helvetica（SCI期刊常用）
plt.rcParams["font.size"] = 12  # 调整字体大小

# 找到最新的OrthoFinder结果目录
results_dir = max(glob.glob(os.path.expanduser("~/workspace/01.mito/02.gb/fasta_output/OrthoFinder/Results_*")), key=os.path.getctime)
gene_count_file = f"{results_dir}/Orthogroups/Orthogroups.GeneCount.tsv"

# 读取基因计数文件
df = pd.read_csv(gene_count_file, sep="\t")

# 计算每次添加株系时的新基因簇数
strains = df.columns[1:-1]
n_strains = len(strains)
n_perms = 100  # 随机排列次数
pangenome_sizes = []

for _ in range(n_perms):
    perm_strains = list(strains)
    shuffle(perm_strains)  # 随机打乱株系顺序
    gene_set = set()
    sizes = []
    for i in range(n_strains):
        new_genes = set(df[df[perm_strains[i]] > 0]["Orthogroup"])
        gene_set.update(new_genes)
        sizes.append(len(gene_set))
    pangenome_sizes.append(sizes)

# 计算平均曲线和标准差（用于误差棒）
pangenome_sizes = np.array(pangenome_sizes)
mean_sizes = pangenome_sizes.mean(axis=0)
std_sizes = pangenome_sizes.std(axis=0)  # 标准差
# 如果要95%置信区间，可以用：1.96 * std_sizes / np.sqrt(n_perms)

# 绘制累积曲线（带误差棒）
plt.figure(figsize=(10, 6))  # 调整图形大小
plt.plot(
    range(1, n_strains + 1),
    mean_sizes,
    marker='o',
    color='#1f77b4',
    linestyle='-',
    linewidth=2,
    markersize=8,
    label="Mean Pangenome Size"
)
plt.fill_between(  # 填充误差范围（标准差）
    range(1, n_strains + 1),
    mean_sizes - std_sizes,
    mean_sizes + std_sizes,
    color='#1f77b4',
    alpha=0.2,
    label="± SD"
)

# 优化图形
plt.xlabel("Number of Strains", fontweight='bold')
plt.ylabel("Number of Orthogroups", fontweight='bold')
plt.title("Pangenome Accumulation Curve for Lentinula edodes mtDNA", fontsize=14, pad=20)
plt.xticks(range(1, n_strains + 1, 1))  # 确保x轴每个株系都有刻度
plt.grid(True, linestyle='--', alpha=0.6)  # 虚线网格
plt.legend(frameon=True, shadow=True)  # 带阴影的图例
plt.tight_layout()  # 防止标签被截断

# 保存高分辨率图片（300 DPI）
output_path = os.path.expanduser("~/workspace/01.mito/02.gb/pangenome_curve.png")
plt.savefig(output_path, dpi=300, bbox_inches='tight')
plt.show()
