#!/usr/bin/env python3
"""
白玉菇线粒体泛基因组变异分析 - 水平一行组合图
所有柱子视觉宽度一致
"""

import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import numpy as np

# 配置
plt.rcParams['font.family'] = ['DejaVu Sans']
plt.rcParams['axes.unicode_minus'] = False
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42

# 字体大小
TITLE_SIZE = 14
LABEL_SIZE = 12
TICK_SIZE = 10
ANNOTATION_SIZE = 10

# 统一柱宽
BAR_WIDTH = 0.65

def create_horizontal_combined():
    """创建水平一行组合图"""
    
    # 柱子数量: 图A=4, 图B=4, 图C=12  
    # 宽度比例 = 4:4:12 = 1:1:3
    fig = plt.figure(figsize=(20, 5))
    gs = gridspec.GridSpec(1, 3, figure=fig, width_ratios=[4, 4, 12], wspace=0.25)
    
    # =========================================================================
    # 图A：变异类型分布
    # =========================================================================
    ax1 = fig.add_subplot(gs[0, 0])
    
    types1 = ['SNPs', 'InDels', 'MNPs', 'Others']
    counts1 = [1563, 522, 246, 175]
    colors1 = ['#5B9BD5', '#C490D1', '#F4B183', '#FF8B8B']
    x1 = np.arange(len(types1))
    
    bars1 = ax1.bar(x1, counts1, width=BAR_WIDTH, color=colors1, edgecolor='black', linewidth=1)
    ax1.set_ylabel('Number of Variants', fontsize=LABEL_SIZE, fontweight='bold')
    ax1.set_xlabel('Variant Type', fontsize=LABEL_SIZE, fontweight='bold')
    ax1.set_title('A. Distribution of Variant Types', fontsize=TITLE_SIZE, fontweight='bold', pad=10, loc='left')
    ax1.set_ylim(0, 1850)
    ax1.set_xlim(-0.5, 3.5)
    ax1.set_xticks(x1)
    ax1.set_xticklabels(types1, fontsize=TICK_SIZE)
    ax1.tick_params(axis='y', labelsize=TICK_SIZE)
    ax1.spines['top'].set_visible(False)
    ax1.spines['right'].set_visible(False)
    
    for bar, count in zip(bars1, counts1):
        ax1.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 30,
                f'{count}', ha='center', va='bottom', fontsize=ANNOTATION_SIZE, fontweight='bold')
    
    # =========================================================================
    # 图B：基因组区域变异分布
    # =========================================================================
    ax2 = fig.add_subplot(gs[0, 1])
    
    regions = ['Intergenic', 'CDS', 'rRNA', 'tRNA']
    counts2 = [1833, 363, 125, 2]
    percentages = ['78.9%', '15.6%', '5.4%', '']
    colors2 = ['#5B9BD5', '#C490D1', '#E6B422', '#DDA0DD']
    x2 = np.arange(len(regions))
    
    bars2 = ax2.bar(x2, counts2, width=BAR_WIDTH, color=colors2, edgecolor='black', linewidth=1)
    ax2.set_ylabel('Number of Variants', fontsize=LABEL_SIZE, fontweight='bold')
    ax2.set_xlabel('Genomic Region', fontsize=LABEL_SIZE, fontweight='bold')
    ax2.set_title('B. Variant Count by Genomic Region', fontsize=TITLE_SIZE, fontweight='bold', pad=10, loc='left')
    ax2.set_ylim(0, 2200)
    ax2.set_xlim(-0.5, 3.5)
    ax2.set_xticks(x2)
    ax2.set_xticklabels(regions, fontsize=TICK_SIZE)
    ax2.tick_params(axis='y', labelsize=TICK_SIZE)
    ax2.spines['top'].set_visible(False)
    ax2.spines['right'].set_visible(False)
    
    for bar, count, pct in zip(bars2, counts2, percentages):
        ax2.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 40,
                f'{count:,}', ha='center', va='bottom', fontsize=ANNOTATION_SIZE, fontweight='bold')
        if pct:
            ax2.text(bar.get_x() + bar.get_width()/2, bar.get_height()/2,
                    pct, ha='center', va='center', fontsize=ANNOTATION_SIZE-1, fontweight='bold', color='white')
    
    # =========================================================================
    # 图C：突变谱
    # =========================================================================
    ax3 = fig.add_subplot(gs[0, 2])
    
    substitutions = ['A>C', 'A>G', 'A>T', 'C>A', 'C>G', 'C>T', 'G>A', 'G>C', 'G>T', 'T>A', 'T>C', 'T>G']
    counts3 = [109, 187, 127, 98, 14, 237, 198, 17, 95, 130, 197, 101]
    
    colors3 = ['#E8875A', '#E8875A', '#E8875A',
               '#5B9BD5', '#5B9BD5', '#5B9BD5',
               '#7FBF7F', '#7FBF7F', '#7FBF7F',
               '#E6B422', '#E6B422', '#E6B422']
    
    x3 = np.arange(len(substitutions))
    bars3 = ax3.bar(x3, counts3, width=BAR_WIDTH, color=colors3, edgecolor='black', linewidth=1)
    
    ax3.set_ylabel('Count', fontsize=LABEL_SIZE, fontweight='bold')
    ax3.set_xlabel('Substitution Type', fontsize=LABEL_SIZE, fontweight='bold')
    ax3.set_title('C. Mutation Spectrum (Ts/Tv = 1.19)', fontsize=TITLE_SIZE, fontweight='bold', pad=10, loc='left')
    ax3.set_ylim(0, 280)
    ax3.set_xlim(-0.5, 11.5)
    ax3.set_xticks(x3)
    ax3.set_xticklabels(substitutions, fontsize=TICK_SIZE)
    ax3.tick_params(axis='y', labelsize=TICK_SIZE)
    ax3.spines['top'].set_visible(False)
    ax3.spines['right'].set_visible(False)
    ax3.grid(axis='y', linestyle='--', alpha=0.3)
    
    for bar, count in zip(bars3, counts3):
        ax3.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 5,
                f'{count}', ha='center', va='bottom', fontsize=ANNOTATION_SIZE-2, fontweight='bold')
    
    # 分隔线
    for i in [2.5, 5.5, 8.5]:
        ax3.axvline(x=i, color='gray', linestyle=':', linewidth=1.5, alpha=0.5)
    
    # 图例
    from matplotlib.patches import Patch
    legend_elements = [
        Patch(facecolor='#E8875A', edgecolor='black', label='A>N'),
        Patch(facecolor='#5B9BD5', edgecolor='black', label='C>N'),
        Patch(facecolor='#7FBF7F', edgecolor='black', label='G>N'),
        Patch(facecolor='#E6B422', edgecolor='black', label='T>N')
    ]
    ax3.legend(handles=legend_elements, loc='upper right', fontsize=9, ncol=4, 
               framealpha=0.9, columnspacing=0.5)
    
    # 保存
    plt.savefig('pangenome_variants_horizontal.pdf', format='pdf', dpi=300, bbox_inches='tight')
    plt.savefig('pangenome_variants_horizontal.png', format='png', dpi=300, bbox_inches='tight')
    plt.close()
    print("✓ 水平组合图已保存: pangenome_variants_horizontal.pdf")

if __name__ == "__main__":
    print("=" * 50)
    print("生成水平一行组合图")
    print("=" * 50)
    create_horizontal_combined()
    print("完成！")
