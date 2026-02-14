import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import glob
import os
from random import shuffle

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

# 计算平均曲线
pangenome_sizes = np.array(pangenome_sizes)
mean_sizes = pangenome_sizes.mean(axis=0)

# 绘制累积曲线
plt.plot(range(1, n_strains + 1), mean_sizes, marker='o', color='#1f77b4')
plt.xlabel("Number of Strains")
plt.ylabel("Number of Orthogroups")
plt.title("Pangenome Accumulation Curve for Lentinula edodes mtDNA")
plt.grid(True)
plt.savefig(os.path.expanduser("~/workspace/01.mito/02.gb/pangenome_curve.png"))
plt.show()
