import os
from Bio import SeqIO
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

# 设置输入目录
input_dir = os.path.expanduser("~/workspace/01.mito/02.gb/1.zuizhong.8.19")

# 目标基因列表
target_genes = ["cox1", "nad5", "cox2", "cob"]

# 统计每个株系中各基因的内含子数量
intron_data = {}
for gb_file in os.listdir(input_dir):
    if gb_file.endswith('.gb'):
        file_path = os.path.join(input_dir, gb_file)
        gene_introns = {gene: 0 for gene in target_genes}  # 初始化计数器
        
        for record in SeqIO.parse(file_path, "genbank"):
            for feature in record.features:
                if feature.type == "intron":
                    # 获取内含子所属的基因
                    intron_gene = feature.qualifiers.get('gene', ['unknown'])[0]
                    # 如果是目标基因，则计数
                    if intron_gene in target_genes:
                        gene_introns[intron_gene] += 1
        
        # 记录数据
        intron_data[gb_file] = gene_introns

# 转换为DataFrame
df = pd.DataFrame.from_dict(intron_data, orient='index')
df = df[target_genes]  # 确保列顺序一致

# 绘制堆叠柱状图
plt.figure(figsize=(12, 6))
df.plot(kind='bar', stacked=True, color=['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728'])

plt.title("Intron Counts by Gene (cox1, nad5, cox2, cob)")
plt.xlabel("Strain")
plt.ylabel("Number of Introns")
plt.xticks(rotation=45)
plt.legend(title="Gene")
plt.tight_layout()

# 保存图片
output_path = os.path.join(input_dir, "intron_counts_by_gene.png")
plt.savefig(output_path, dpi=300, bbox_inches='tight')
plt.show()
