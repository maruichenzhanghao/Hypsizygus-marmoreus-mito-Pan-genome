import os
from Bio import SeqIO
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import gaussian_kde

# 获取当前目录下的所有 .fasta 文件
fasta_files = [f for f in os.listdir('.') if f.endswith('.fasta')]

# 存储每个 .fasta 文件的序列长度
sequence_lengths = []

for fasta_file in fasta_files:
    for record in SeqIO.parse(fasta_file, "fasta"):
        sequence_lengths.append(len(record.seq))

# 将序列长度转换为 kb 单位
sequence_lengths_kb = [length / 1000 for length in sequence_lengths]

# 绘制直方图和密度曲线
fig, ax1 = plt.subplots()

# 直方图
bins = np.arange(min(sequence_lengths_kb), max(sequence_lengths_kb) + 2, 1)  # 根据你的数据范围调整 bin 的范围和间隔
ax1.hist(sequence_lengths_kb, bins=bins, color='gray', edgecolor='black', alpha=0.7)
ax1.set_xlabel('Mitogenome size (kb)')
ax1.set_ylabel('Count')
ax1.set_title('Mitogenome Size Distribution')

# 密度曲线
density = gaussian_kde(sequence_lengths_kb)
x_range = np.linspace(min(sequence_lengths_kb), max(sequence_lengths_kb), 100)
ax2 = ax1.twinx()
ax2.plot(x_range, density(x_range), 'k--')
ax2.set_ylabel('Density')

# 显示网格
plt.grid(True)

# 显示图形
plt.show()
