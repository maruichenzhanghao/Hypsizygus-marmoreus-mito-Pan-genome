#!/usr/bin/env python3
"""
白玉菇线粒体泛基因组变异分析图表
基于minigraph-cactus分析结果
输出矢量PDF格式
"""

import matplotlib.pyplot as plt
import matplotlib.font_manager as fm
import numpy as np
import os

# 配置中文字体
font_path = '/usr/share/fonts/opentype/noto/NotoSansCJK-Regular.ttc'
if os.path.exists(font_path):
    fm.fontManager.addfont(font_path)
    plt.rcParams['font.family'] = ['Noto Sans CJK SC', 'DejaVu Sans']
else:
    plt.rcParams['font.family'] = ['DejaVu Sans']

plt.rcParams['axes.unicode_minus'] = False
plt.rcParams['pdf.fonttype'] = 42  # 确保PDF中字体可编辑
plt.rcParams['ps.fonttype'] = 42

# =============================================================================
# 图1: 变异类型分布 + Ts/Tv比值
# =============================================================================
def plot_variant_types_and_tstv():
    """绑制变异类型分布和Ts/Tv比值图"""
    # 使用gridspec控制子图宽度比例，使柱子视觉宽度一致
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5), 
                                    gridspec_kw={'width_ratios': [4, 2]})
    
    # 统一柱宽（相对于各自x轴范围）
    bar_width = 0.6
    
    # 左图：变异类型分布（4个柱子）
    types = ['SNPs', 'InDels', 'MNPs', 'Others']
    counts = [1563, 522, 246, 175]
    colors = ['#5B9BD5', '#C490D1', '#F4B183', '#FF8B8B']
    x1 = np.arange(len(types))
    
    bars1 = ax1.bar(x1, counts, width=bar_width, color=colors, edgecolor='black', linewidth=0.8)
    ax1.set_ylabel('Number of Variants', fontsize=12, fontweight='bold')
    ax1.set_xlabel('Type', fontsize=12, fontweight='bold')
    ax1.set_title('Distribution of Variant Types', fontsize=14, fontweight='bold', pad=15)
    ax1.set_ylim(0, 1800)
    ax1.set_xlim(-0.5, 3.5)
    ax1.set_xticks(x1)
    ax1.set_xticklabels(types)
    ax1.spines['top'].set_visible(False)
    ax1.spines['right'].set_visible(False)
    
    # 添加数值标签
    for bar, count in zip(bars1, counts):
        ax1.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 30,
                f'{count}', ha='center', va='bottom', fontsize=11, fontweight='bold')
    
    # 右图：Ts/Tv比值（2个柱子）
    tstv_types = ['Transitions', 'Transversions']
    tstv_counts = [820, 690]
    tstv_colors = ['#5B9BD5', '#F08080']
    x2 = np.arange(len(tstv_types))
    
    bars2 = ax2.bar(x2, tstv_counts, width=bar_width, color=tstv_colors, edgecolor='black', linewidth=0.8)
    ax2.set_ylabel('Count', fontsize=12, fontweight='bold')
    ax2.set_xlabel('Type', fontsize=12, fontweight='bold')
    ax2.set_title('Ts/Tv Ratio = 1.19', fontsize=14, fontweight='bold', pad=15)
    ax2.set_ylim(0, 1000)
    ax2.set_xlim(-0.5, 1.5)
    ax2.set_xticks(x2)
    ax2.set_xticklabels(tstv_types)
    ax2.spines['top'].set_visible(False)
    ax2.spines['right'].set_visible(False)
    
    # 添加数值标签
    for bar, count in zip(bars2, tstv_counts):
        ax2.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 20,
                f'{count}', ha='center', va='bottom', fontsize=11, fontweight='bold')
    
    plt.tight_layout()
    plt.savefig('fig1_variant_types_tstv.pdf', format='pdf', dpi=300, bbox_inches='tight')
    plt.savefig('fig1_variant_types_tstv.png', format='png', dpi=300, bbox_inches='tight')
    plt.close()
    print("✓ 图1已保存: fig1_variant_types_tstv.pdf")

# =============================================================================
# 图2: 基因组区域变异分布 + 变异密度
# =============================================================================
def plot_genomic_region_variants():
    """绑制基因组区域变异分布图"""
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13, 5.5))
    
    # 统一柱宽
    bar_width = 0.6
    
    # 数据
    regions = ['Intergenic', 'CDS', 'rRNA', 'tRNA']
    counts = [1833, 363, 125, 2]
    percentages = ['78.9%', '15.6%', '5.4%', '']
    densities = [36.90, 6.88, 7.51, 0.99]
    enrichments = ['1.9×', '0.4×', '0.4×', '']
    colors = ['#5B9BD5', '#C490D1', '#E6B422', '#DDA0DD']
    x = np.arange(len(regions))
    
    # 左图：变异计数
    bars1 = ax1.bar(x, counts, width=bar_width, color=colors, edgecolor='black', linewidth=0.8)
    ax1.set_ylabel('Number of Variants', fontsize=12, fontweight='bold')
    ax1.set_xlabel('Genomic Region', fontsize=12, fontweight='bold')
    ax1.set_title('Variant Count by Genomic Region', fontsize=14, fontweight='bold', pad=15)
    ax1.set_ylim(0, 2200)
    ax1.set_xticks(x)
    ax1.set_xticklabels(regions)
    ax1.spines['top'].set_visible(False)
    ax1.spines['right'].set_visible(False)
    
    # 添加总数标注框
    ax1.text(0.95, 0.95, 'Total: 2,323', transform=ax1.transAxes, fontsize=10,
             verticalalignment='top', horizontalalignment='right',
             bbox=dict(boxstyle='round', facecolor='lightgray', alpha=0.8))
    
    # 添加数值和百分比标签
    for bar, count, pct in zip(bars1, counts, percentages):
        ax1.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 30,
                f'{count:,}', ha='center', va='bottom', fontsize=11, fontweight='bold')
        if pct:
            ax1.text(bar.get_x() + bar.get_width()/2, bar.get_height()/2,
                    pct, ha='center', va='center', fontsize=10, fontweight='bold', color='white')
    
    # 右图：变异密度
    bars2 = ax2.bar(x, densities, width=bar_width, color=colors, edgecolor='black', linewidth=0.8)
    ax2.set_ylabel('Variant Density (per kb)', fontsize=12, fontweight='bold')
    ax2.set_xlabel('Genomic Region', fontsize=12, fontweight='bold')
    ax2.set_title('Variant Density by Genomic Region', fontsize=14, fontweight='bold', pad=15)
    ax2.set_ylim(0, 45)
    ax2.set_xticks(x)
    ax2.set_xticklabels(regions)
    ax2.spines['top'].set_visible(False)
    ax2.spines['right'].set_visible(False)
    
    # 添加平均值标注框
    ax2.text(0.95, 0.95, 'Avg: 13.07 per kb', transform=ax2.transAxes, fontsize=10,
             verticalalignment='top', horizontalalignment='right',
             bbox=dict(boxstyle='round', facecolor='lightyellow', alpha=0.8))
    
    # 添加副标题
    ax2.text(0.5, 1.02, 'Numbers = Enrichment factor', transform=ax2.transAxes, fontsize=9,
             ha='center', color='gray', style='italic')
    
    # 添加数值和富集因子标签
    for bar, density, enrich in zip(bars2, densities, enrichments):
        ax2.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.8,
                f'{density:.2f}', ha='center', va='bottom', fontsize=11, fontweight='bold')
        if enrich:
            ax2.text(bar.get_x() + bar.get_width()/2, bar.get_height()/2,
                    enrich, ha='center', va='center', fontsize=10, fontweight='bold', color='white')
    
    plt.tight_layout()
    plt.savefig('fig2_genomic_region_variants.pdf', format='pdf', dpi=300, bbox_inches='tight')
    plt.savefig('fig2_genomic_region_variants.png', format='png', dpi=300, bbox_inches='tight')
    plt.close()
    print("✓ 图2已保存: fig2_genomic_region_variants.pdf")

# =============================================================================
# 图3: 突变谱 - 碱基替换分布
# =============================================================================
def plot_mutation_spectrum():
    """绘制突变谱图"""
    fig, ax = plt.subplots(figsize=(14, 6))
    
    # 统一柱宽
    bar_width = 0.6
    
    # 数据 - 按碱基替换类型排列
    substitutions = ['A>C', 'A>G', 'A>T', 'C>A', 'C>G', 'C>T', 'G>A', 'G>C', 'G>T', 'T>A', 'T>C', 'T>G']
    counts = [109, 187, 127, 98, 14, 237, 198, 17, 95, 130, 197, 101]
    percentages = ['7.2%', '12.4%', '8.4%', '6.5%', '', '15.7%', '13.1%', '', '6.3%', '8.6%', '13.0%', '6.7%']
    
    # 按来源碱基分组设置颜色
    # A>N: 橙色, C>N: 蓝色, G>N: 绿色, T>N: 黄色
    colors = ['#E8875A', '#E8875A', '#E8875A',  # A>N
              '#5B9BD5', '#5B9BD5', '#5B9BD5',  # C>N
              '#7FBF7F', '#7FBF7F', '#7FBF7F',  # G>N
              '#E6B422', '#E6B422', '#E6B422']  # T>N
    
    # 绘制柱状图
    x = np.arange(len(substitutions))
    bars = ax.bar(x, counts, color=colors, edgecolor='black', linewidth=0.8, width=bar_width)
    
    ax.set_ylabel('Count', fontsize=12, fontweight='bold')
    ax.set_xlabel('Substitution Type', fontsize=12, fontweight='bold')
    ax.set_title('Mutation Spectrum: Distribution of Base Substitutions', fontsize=14, fontweight='bold', pad=15)
    ax.set_xticks(x)
    ax.set_xticklabels(substitutions, fontsize=10)
    ax.set_ylim(0, 280)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.grid(axis='y', linestyle='--', alpha=0.3)
    
    # 添加数值和百分比标签
    for bar, count, pct in zip(bars, counts, percentages):
        ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 5,
                f'{count}', ha='center', va='bottom', fontsize=10, fontweight='bold')
        if pct:
            ax.text(bar.get_x() + bar.get_width()/2, bar.get_height()/2,
                    pct, ha='center', va='center', fontsize=8, fontweight='bold',
                    color='white', bbox=dict(boxstyle='round,pad=0.2', facecolor='gray', alpha=0.6))
    
    # 添加图例
    from matplotlib.patches import Patch
    legend_elements = [
        Patch(facecolor='#E8875A', edgecolor='black', label='A>N'),
        Patch(facecolor='#5B9BD5', edgecolor='black', label='C>N'),
        Patch(facecolor='#7FBF7F', edgecolor='black', label='G>N'),
        Patch(facecolor='#E6B422', edgecolor='black', label='T>N')
    ]
    ax.legend(handles=legend_elements, title='Mutation Category', loc='upper left', 
              framealpha=0.9, fontsize=9, title_fontsize=10)
    
    # 添加统计信息标注框
    stats_text = 'Total: 1,510'
    ax.text(0.98, 0.95, stats_text, transform=ax.transAxes, fontsize=10,
            verticalalignment='top', horizontalalignment='right',
            bbox=dict(boxstyle='round', facecolor='lightgray', alpha=0.8))
    
    # Ts/Tv信息框
    tstv_text = 'Ts/Tv: 1.19\nTs: 819 | Tv: 691'
    ax.text(0.98, 0.82, tstv_text, transform=ax.transAxes, fontsize=9,
            verticalalignment='top', horizontalalignment='right',
            bbox=dict(boxstyle='round', facecolor='lightgreen', alpha=0.7))
    
    # 最频繁替换标注
    ax.text(0.98, 0.65, 'Most frequent:\nC>T (237)', transform=ax.transAxes, fontsize=9,
            verticalalignment='top', horizontalalignment='right',
            bbox=dict(boxstyle='round', facecolor='lightyellow', alpha=0.8))
    
    # 添加垂直分隔线
    for i in [2.5, 5.5, 8.5]:
        ax.axvline(x=i, color='gray', linestyle=':', linewidth=1, alpha=0.5)
    
    plt.tight_layout()
    plt.savefig('fig3_mutation_spectrum.pdf', format='pdf', dpi=300, bbox_inches='tight')
    plt.savefig('fig3_mutation_spectrum.png', format='png', dpi=300, bbox_inches='tight')
    plt.close()
    print("✓ 图3已保存: fig3_mutation_spectrum.pdf")

# =============================================================================
# 主函数
# =============================================================================
def main():
    print("=" * 50)
    print("白玉菇线粒体泛基因组变异分析图表生成")
    print("=" * 50)
    print()
    
    # 绑制所有图表
    plot_variant_types_and_tstv()
    plot_genomic_region_variants()
    plot_mutation_spectrum()
    
    print()
    print("=" * 50)
    print("所有图表已生成完毕！")
    print("PDF文件（矢量格式）：")
    print("  - fig1_variant_types_tstv.pdf")
    print("  - fig2_genomic_region_variants.pdf")
    print("  - fig3_mutation_spectrum.pdf")
    print("=" * 50)

if __name__ == "__main__":
    main()
