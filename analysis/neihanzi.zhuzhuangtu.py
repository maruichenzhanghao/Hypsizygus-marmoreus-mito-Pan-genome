import os
from Bio import SeqIO
import matplotlib.pyplot as plt
import pandas as pd

# 设置输入目录
input_dir = os.path.expanduser("~/workspace/01.mito/02.gb/1.gb")

# 统计内含子数量和类型
intron_data = {}
for gb_file in os.listdir(input_dir):
    if gb_file.endswith('.gb'):  # 检查文件扩展名
        file_path = os.path.join(input_dir, gb_file)  # 正确拼接路径
        for record in SeqIO.parse(file_path, "genbank"):
            intron_count = 0
            intron_types = {}
            for feature in record.features:
                if feature.type == "intron":
                    intron_count += 1
                    # 尝试提取内含子类型（如 Group I, Group IB 等）
                    intron_type = feature.qualifiers.get('note', ['unknown'])[0]
                    intron_types[intron_type] = intron_types.get(intron_type, 0) + 1
            intron_data[gb_file] = {
                'total_introns': intron_count,
                'intron_types': intron_types
            }

# 转换为 DataFrame 便于分析
df = pd.DataFrame.from_dict(intron_data, orient='index')

# 1. 绘制总内含子数量的柱状图
plt.figure(figsize=(10, 6))
df['total_introns'].plot(kind='bar', color='#1f77b4')
plt.title("Total Intron Counts per Strain")
plt.xlabel("Strain")
plt.ylabel("Number of Introns")
plt.xticks(rotation=45)
plt.tight_layout()
plt.savefig("total_intron_counts.png", dpi=300)
plt.show()

# 2. 绘制内含子类型分布（以第一个株系为例）
if not df.empty:
    example_strain = df.index[0]
    intron_types = df.loc[example_strain, 'intron_types']
    
    plt.figure(figsize=(10, 6))
    pd.Series(intron_types).plot(kind='bar', color='#ff7f0e')
    plt.title(f"Intron Types in {example_strain}")
    plt.xlabel("Intron Type")
    plt.ylabel("Count")
    plt.xticks(rotation=45)
    plt.tight_layout()
    plt.savefig(f"intron_types_{example_strain}.png", dpi=300)
    plt.show()
