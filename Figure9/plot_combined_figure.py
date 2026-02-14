#!/usr/bin/env python3
"""
白玉菇线粒体泛基因组变异分析 - 组合图
三个图表左右排列，统一字体和柱宽
"""

import matplotlib.pyplot as plt
import matplotlib.font_manager as fm
import matplotlib.gridspec as gridspec
import numpy as np
import os

# 配置字体
plt.rcParams['font.family'] = ['DejaVu Sans']
plt.rcParams['axes.unicode_minus'] = False
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42

# 全局字体大小设置（放大）
TITLE_SIZE = 16
LABEL_SIZE = 14
TICK_SIZE = 12
ANNOTATION_SIZE = 12
LEGEND_SIZE = 11

# 统一柱宽
BAR_WIDTH = 0.55

def create_combined_figure():
    """创建组合图"""
    
    # 创建3行2列的布局
    fig = plt.figure(figsize=(16, 18))
    
    # 使用GridSpec精确控制子图大小
    # 每行左右两个子图的宽度比例根据柱子数量调整
    gs = gridspec.GridSpec(3, 2, figure=fig, 
                           width_ratios=[2, 1],  # 左:右 = 2:1
                           height_ratios=[1, 1, 1.2],  # 第三行稍高
                           hspace=0.35, wspace=0.25)
    
    # =========================================================================
    # 第一行：变异类型分布 + Ts/Tv比值
    # =========================================================================
    ax1 = fig.add_subplot(gs[0, 0])
    ax2 = fig.add_subplot(gs[0, 1])
    
    # 左图：变异类型分布（4个柱子）
    types1 = ['SNPs', 'InDels', 'MNPs', 'Others']
    counts1 = [1563, 522, 246, 175]
    colors1 = ['#5B9BD5', '#C490D1', '#F4B183', '#FF8B8B']
    x1 = np.arange(len(types1))
    
    bars1 = ax1.bar(x1, counts1, width=BAR_WIDTH, color=colors1, edgecolor='black', linewidth=1)
    ax1.set_ylabel('Number of Variants', fontsize=LABEL_SIZE, fontweight='bold')
    ax1.set_xlabel('Type', fontsize=LABEL_SIZE, fontweight='bold')
    ax1.set_title('A. Distribution of Variant Types', fontsize=TITLE_SIZE, fontweight='bold', pad=15, loc='left')
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
    
    # 右图：Ts/Tv比值（2个柱子）
    types2 = ['Transitions', 'Transversions']
    counts2 = [820, 690]
    colors2 = ['#5B9BD5', '#F08080']
    x2 = np.arange(len(types2))
    
    bars2 = ax2.bar(x2, counts2, width=BAR_WIDTH, color=colors2, edgecolor='black', linewidth=1)
    ax2.set_ylabel('Count', fontsize=LABEL_SIZE, fontweight='bold')
    ax2.set_xlabel('Type', fontsize=LABEL_SIZE, fontweight='bold')
    ax2.set_title('B. Ts/Tv Ratio = 1.19', fontsize=TITLE_SIZE, fontweight='bold', pad=15, loc='left')
    ax2.set_ylim(0, 1000)
    ax2.set_xlim(-0.5, 1.5)
    ax2.set_xticks(x2)
    ax2.set_xticklabels(types2, fontsize=TICK_SIZE)
    ax2.tick_params(axis='y', labelsize=TICK_SIZE)
    ax2.spines['top'].set_visible(False)
    ax2.spines['right'].set_visible(False)
    
    for bar, count in zip(bars2, counts2):
        ax2.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 20,
                f'{count}', ha='center', va='bottom', fontsize=ANNOTATION_SIZE, fontweight='bold')
    
    # =========================================================================
    # 第二行：基因组区域变异分布
    # =========================================================================
    ax3 = fig.add_subplot(gs[1, 0])
    ax4 = fig.add_subplot(gs[1, 1])
    
    # 数据
    regions = ['Intergenic', 'CDS', 'rRNA', 'tRNA']
    counts3 = [1833, 363, 125, 2]
    percentages = ['78.9%', '15.6%', '5.4%', '']
    densities = [36.90, 6.88, 7.51, 0.99]
    enrichments = ['1.9×', '0.4×', '0.4×', '']
    colors3 = ['#5B9BD5', '#C490D1', '#E6B422', '#DDA0DD']
    x3 = np.arange(len(regions))
    
    # 左图：变异计数（4个柱子）
    bars3 = ax3.bar(x3, counts3, width=BAR_WIDTH, color=colors3, edgecolor='black', linewidth=1)
    ax3.set_ylabel('Number of Variants', fontsize=LABEL_SIZE, fontweight='bold')
    ax3.set_xlabel('Genomic Region', fontsize=LABEL_SIZE, fontweight='bold')
    ax3.set_title('C. Variant Count by Genomic Region', fontsize=TITLE_SIZE, fontweight='bold', pad=15, loc='left')
    ax3.set_ylim(0, 2200)
    ax3.set_xlim(-0.5, 3.5)
    ax3.set_xticks(x3)
    ax3.set_xticklabels(regions, fontsize=TICK_SIZE)
    ax3.tick_params(axis='y', labelsize=TICK_SIZE)
    ax3.spines['top'].set_visible(False)
    ax3.spines['right'].set_visible(False)
    
    ax3.text(0.95, 0.95, 'Total: 2,323', transform=ax3.transAxes, fontsize=LEGEND_SIZE,
             verticalalignment='top', horizontalalignment='right',
             bbox=dict(boxstyle='round', facecolor='lightgray', alpha=0.8))
    
    for bar, count, pct in zip(bars3, counts3, percentages):
        ax3.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 40,
                f'{count:,}', ha='center', va='bottom', fontsize=ANNOTATION_SIZE, fontweight='bold')
        if pct:
            ax3.text(bar.get_x() + bar.get_width()/2, bar.get_height()/2,
                    pct, ha='center', va='center', fontsize=ANNOTATION_SIZE-1, fontweight='bold', color='white')
    
    # 右图：变异密度（4个柱子，但用2:1比例的右子图）
    # 为了保持柱宽一致，这里也使用4个数据点
    ax4_2 = fig.add_subplot(gs[1, 1])
    
    # 只显示前两个密度较高的区域，便于视觉对比
    bars4 = ax4_2.bar(x3, densities, width=BAR_WIDTH, color=colors3, edgecolor='black', linewidth=1)
    ax4_2.set_ylabel('Variant Density (per kb)', fontsize=LABEL_SIZE, fontweight='bold')
    ax4_2.set_xlabel('Genomic Region', fontsize=LABEL_SIZE, fontweight='bold')
    ax4_2.set_title('D. Variant Density by Region', fontsize=TITLE_SIZE, fontweight='bold', pad=15, loc='left')
    ax4_2.set_ylim(0, 45)
    ax4_2.set_xlim(-0.5, 3.5)
    ax4_2.set_xticks(x3)
    ax4_2.set_xticklabels(regions, fontsize=TICK_SIZE-1, rotation=15, ha='right')
    ax4_2.tick_params(axis='y', labelsize=TICK_SIZE)
    ax4_2.spines['top'].set_visible(False)
    ax4_2.spines['right'].set_visible(False)
    
    ax4_2.text(0.95, 0.95, 'Avg: 13.07/kb', transform=ax4_2.transAxes, fontsize=LEGEND_SIZE,
             verticalalignment='top', horizontalalignment='right',
             bbox=dict(boxstyle='round', facecolor='lightyellow', alpha=0.8))
    
    for bar, density, enrich in zip(bars4, densities, enrichments):
        ax4_2.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 1,
                f'{density:.1f}', ha='center', va='bottom', fontsize=ANNOTATION_SIZE-1, fontweight='bold')
    
    # =========================================================================
    # 第三行：突变谱（横跨整行）
    # =========================================================================
    ax5 = fig.add_subplot(gs[2, :])
    
    substitutions = ['A>C', 'A>G', 'A>T', 'C>A', 'C>G', 'C>T', 'G>A', 'G>C', 'G>T', 'T>A', 'T>C', 'T>G']
    counts5 = [109, 187, 127, 98, 14, 237, 198, 17, 95, 130, 197, 101]
    percentages5 = ['7.2%', '12.4%', '8.4%', '6.5%', '', '15.7%', '13.1%', '', '6.3%', '8.6%', '13.0%', '6.7%']
    
    colors5 = ['#E8875A', '#E8875A', '#E8875A',
               '#5B9BD5', '#5B9BD5', '#5B9BD5',
               '#7FBF7F', '#7FBF7F', '#7FBF7F',
               '#E6B422', '#E6B422', '#E6B422']
    
    x5 = np.arange(len(substitutions))
    bars5 = ax5.bar(x5, counts5, color=colors5, edgecolor='black', linewidth=1, width=BAR_WIDTH)
    
    ax5.set_ylabel('Count', fontsize=LABEL_SIZE, fontweight='bold')
    ax5.set_xlabel('Substitution Type', fontsize=LABEL_SIZE, fontweight='bold')
    ax5.set_title('E. Mutation Spectrum: Distribution of Base Substitutions', fontsize=TITLE_SIZE, fontweight='bold', pad=15, loc='left')
    ax5.set_xticks(x5)
    ax5.set_xticklabels(substitutions, fontsize=TICK_SIZE)
    ax5.tick_params(axis='y', labelsize=TICK_SIZE)
    ax5.set_ylim(0, 280)
    ax5.spines['top'].set_visible(False)
    ax5.spines['right'].set_visible(False)
    ax5.grid(axis='y', linestyle='--', alpha=0.3)
    
    for bar, count, pct in zip(bars5, counts5, percentages5):
        ax5.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 5,
                f'{count}', ha='center', va='bottom', fontsize=ANNOTATION_SIZE-1, fontweight='bold')
        if pct:
            ax5.text(bar.get_x() + bar.get_width()/2, bar.get_height()/2,
                    pct, ha='center', va='center', fontsize=9, fontweight='bold',
                    color='white', bbox=dict(boxstyle='round,pad=0.15', facecolor='gray', alpha=0.6))
    
    # 图例
    from matplotlib.patches import Patch
    legend_elements = [
        Patch(facecolor='#E8875A', edgecolor='black', label='A>N'),
        Patch(facecolor='#5B9BD5', edgecolor='black', label='C>N'),
        Patch(facecolor='#7FBF7F', edgecolor='black', label='G>N'),
        Patch(facecolor='#E6B422', edgecolor='black', label='T>N')
    ]
    ax5.legend(handles=legend_elements, title='Mutation Category', loc='upper left', 
              framealpha=0.9, fontsize=LEGEND_SIZE, title_fontsize=LEGEND_SIZE)
    
    # 统计信息框
    ax5.text(0.98, 0.95, 'Total: 1,510', transform=ax5.transAxes, fontsize=LEGEND_SIZE,
            verticalalignment='top', horizontalalignment='right',
            bbox=dict(boxstyle='round', facecolor='lightgray', alpha=0.8))
    
    ax5.text(0.98, 0.82, 'Ts/Tv: 1.19\nTs: 819 | Tv: 691', transform=ax5.transAxes, fontsize=LEGEND_SIZE,
            verticalalignment='top', horizontalalignment='right',
            bbox=dict(boxstyle='round', facecolor='lightgreen', alpha=0.7))
    
    # 分隔线
    for i in [2.5, 5.5, 8.5]:
        ax5.axvline(x=i, color='gray', linestyle=':', linewidth=1.5, alpha=0.5)
    
    # 保存
    plt.savefig('pangenome_variants_combined.pdf', format='pdf', dpi=300, bbox_inches='tight')
    plt.savefig('pangenome_variants_combined.png', format='png', dpi=300, bbox_inches='tight')
    plt.close()
    print("✓ 组合图已保存: pangenome_variants_combined.pdf")

if __name__ == "__main__":
    print("=" * 50)
    print("生成泛基因组变异分析组合图")
    print("=" * 50)
    create_combined_figure()
    print("=" * 50)
    print("完成！")
