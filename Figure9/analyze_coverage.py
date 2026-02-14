#!/usr/bin/env python3
"""
线粒体覆盖深度分析：
1. 统计原始测序量（total raw reads & bases）
2. 检测断裂/低覆盖区域
3. 绘制每个样本的覆盖深度折线图
"""
import os
import glob
import subprocess
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
import pandas as pd

WORKDIR = '/home/maxinxin/workspace/01.mito/03.pangenome.cactus/coverage_analysis'
DEPTH_DIR = os.path.join(WORKDIR, 'depth_data')
PLOT_DIR = os.path.join(WORKDIR, 'plots')
HIFI_DIR = '/home/maxinxin/data/02.pangenome/01.zicexu/01.hifi'
ILLUMINA_DIR = '/home/maxinxin/workspace/03.pangenome/01.2dai/new79/02.yiwancheng/completed_fastq'

os.makedirs(PLOT_DIR, exist_ok=True)

# ============================================================
# Sample definitions
# ============================================================
HIFI_SAMPLES = {
    'f2': f'{HIFI_DIR}/F2/F2.hifi_reads.fastq.gz',
    'f4': f'{HIFI_DIR}/F4/F4.hifi_reads.fastq.gz',
    'nn12-1': f'{HIFI_DIR}/NN12-1/NN121.hifi_reads.fastq.gz',
    'nn12-17': f'{HIFI_DIR}/NN12-17/NN1217.hifi_reads.fastq.gz',
}

SRR12151_SAMPLES = ['SRR12151860', 'SRR12151871', 'SRR12151875', 'SRR12151883']
SRR8699_SAMPLES = ['SRR8699796', 'SRR8699797', 'SRR8699800', 'SRR8699801',
                   'SRR8699802', 'SRR8699803', 'SRR8699804', 'SRR8699805',
                   'SRR8699808', 'SRR8699809', 'SRR8699811', 'SRR8699813',
                   'SRR8699814', 'SRR8699815', 'SRR8699816', 'SRR8699817',
                   'SRR8699833', 'SRR8699834', 'SRR8699835', 'SRR8699837']

ILLUMINA_SAMPLES = {s: (f'{ILLUMINA_DIR}/{s}_1.trimmed.fastq.gz',
                         f'{ILLUMINA_DIR}/{s}_2.trimmed.fastq.gz')
                    for s in SRR12151_SAMPLES + SRR8699_SAMPLES}

# ============================================================
# 1. 统计原始测序量
# ============================================================
def count_fastq_stats(fastq_path):
    """使用 seqkit/awk 快速统计 fastq.gz 的 reads 数和 总碱基数"""
    cmd = f"zcat '{fastq_path}' | awk 'NR%4==2{{n++;b+=length($0)}}END{{print n, b}}'"
    result = subprocess.run(cmd, shell=True, capture_output=True, text=True, timeout=600)
    parts = result.stdout.strip().split()
    if len(parts) == 2:
        return int(parts[0]), int(parts[1])
    return 0, 0

print("=" * 60)
print("统计原始测序量...")
print("=" * 60)

raw_stats = {}

# HiFi samples
for sample, fq in HIFI_SAMPLES.items():
    print(f"  {sample} (HiFi)...")
    reads, bases = count_fastq_stats(fq)
    raw_stats[sample] = {'seq_type': 'hifi', 'raw_reads': reads, 'raw_bases': bases}
    print(f"    Raw reads: {reads:,}, Raw bases: {bases:,} ({bases/1e9:.2f} Gb)")

# Illumina samples (count both R1+R2)
for sample, (fq1, fq2) in ILLUMINA_SAMPLES.items():
    print(f"  {sample} (Illumina)...")
    r1, b1 = count_fastq_stats(fq1)
    r2, b2 = count_fastq_stats(fq2)
    total_reads = r1 + r2
    total_bases = b1 + b2
    raw_stats[sample] = {'seq_type': 'illumina', 'raw_reads': total_reads, 'raw_bases': total_bases}
    print(f"    Raw reads: {total_reads:,}, Raw bases: {total_bases:,} ({total_bases/1e9:.2f} Gb)")

# ============================================================
# 2. 读取深度数据 + 检测断裂区域
# ============================================================
print("\n" + "=" * 60)
print("分析覆盖深度 & 检测断裂区域...")
print("=" * 60)

LOW_COV_THRESHOLD = 10  # 低于此深度视为"低覆盖"
ZERO_COV_THRESHOLD = 0  # 0覆盖视为"断裂"

all_depth_data = {}
gap_report = {}

depth_files = sorted(glob.glob(os.path.join(DEPTH_DIR, '*.depth.tsv')))

for df in depth_files:
    sample = os.path.basename(df).replace('.depth.tsv', '')
    
    # Read depth data
    positions = []
    depths = []
    with open(df) as f:
        for line in f:
            parts = line.strip().split('\t')
            positions.append(int(parts[1]))
            depths.append(int(parts[2]))
    
    positions = np.array(positions)
    depths = np.array(depths)
    all_depth_data[sample] = (positions, depths)
    
    genome_len = len(positions)
    
    # Find zero coverage regions (gaps/breaks)
    zero_mask = depths == 0
    zero_count = np.sum(zero_mask)
    
    # Find low coverage regions (<10X)
    low_mask = depths < LOW_COV_THRESHOLD
    low_count = np.sum(low_mask)
    
    # Identify contiguous zero-coverage regions
    zero_regions = []
    if zero_count > 0:
        in_gap = False
        gap_start = 0
        for i in range(len(depths)):
            if depths[i] == 0 and not in_gap:
                in_gap = True
                gap_start = positions[i]
            elif depths[i] > 0 and in_gap:
                in_gap = False
                zero_regions.append((gap_start, positions[i-1], positions[i-1] - gap_start + 1))
        if in_gap:
            zero_regions.append((gap_start, positions[-1], positions[-1] - gap_start + 1))
    
    # Similarly for low coverage
    low_regions = []
    if low_count > 0:
        in_low = False
        low_start = 0
        for i in range(len(depths)):
            if depths[i] < LOW_COV_THRESHOLD and not in_low:
                in_low = True
                low_start = positions[i]
            elif depths[i] >= LOW_COV_THRESHOLD and in_low:
                in_low = False
                low_regions.append((low_start, positions[i-1], positions[i-1] - low_start + 1))
        if in_low:
            low_regions.append((low_start, positions[-1], positions[-1] - low_start + 1))
    
    gap_report[sample] = {
        'genome_len': genome_len,
        'zero_bp': zero_count,
        'zero_regions': zero_regions,
        'low_bp': low_count,
        'low_regions': low_regions,
    }
    
    has_gap = "YES" if zero_count > 0 else "NO"
    print(f"  {sample:20s}: 0X={zero_count}bp ({len(zero_regions)} regions), "
          f"<{LOW_COV_THRESHOLD}X={low_count}bp ({len(low_regions)} regions) → 断裂: {has_gap}")

# ============================================================
# 3. 更新 TSV 文件
# ============================================================
print("\n" + "=" * 60)
print("更新统计表...")
print("=" * 60)

# Read existing TSV
tsv_path = os.path.join(WORKDIR, 'mito_coverage_stats.tsv')
df_existing = pd.read_csv(tsv_path, sep='\t')

# Add new columns
df_existing['Raw_Total_Reads'] = df_existing['Sample'].map(
    lambda s: raw_stats.get(s, {}).get('raw_reads', 'NA'))
df_existing['Raw_Total_Bases'] = df_existing['Sample'].map(
    lambda s: raw_stats.get(s, {}).get('raw_bases', 'NA'))
df_existing['Raw_Total_Bases_Gb'] = df_existing['Sample'].map(
    lambda s: f"{raw_stats.get(s, {}).get('raw_bases', 0)/1e9:.2f}" 
    if s in raw_stats else 'NA')
df_existing['Zero_Coverage_bp'] = df_existing['Sample'].map(
    lambda s: gap_report.get(s, {}).get('zero_bp', 'NA'))
df_existing['Zero_Coverage_Regions'] = df_existing['Sample'].map(
    lambda s: len(gap_report.get(s, {}).get('zero_regions', [])))
df_existing['Has_Breakage'] = df_existing['Sample'].map(
    lambda s: 'YES' if gap_report.get(s, {}).get('zero_bp', 0) > 0 else 'NO')
df_existing['Low_Cov_bp(<10X)'] = df_existing['Sample'].map(
    lambda s: gap_report.get(s, {}).get('low_bp', 'NA'))

# Save updated TSV
df_existing.to_csv(tsv_path, sep='\t', index=False)
print(f"已更新: {tsv_path}")

# ============================================================
# 4. 绘制覆盖深度图
# ============================================================
print("\n" + "=" * 60)
print("绘制覆盖深度图...")
print("=" * 60)

# --- 4a. Individual plots ---
for sample in sorted(all_depth_data.keys()):
    positions, depths = all_depth_data[sample]
    seq_type = raw_stats.get(sample, {}).get('seq_type', 'unknown')
    mean_d = np.mean(depths)
    median_d = np.median(depths)
    
    fig, ax = plt.subplots(figsize=(14, 4))
    
    # Use window average for smoother plot (100bp windows)
    window = 100
    if len(depths) > window:
        n_windows = len(depths) // window
        pos_avg = np.array([positions[i*window] for i in range(n_windows)])
        depth_avg = np.array([np.mean(depths[i*window:(i+1)*window]) for i in range(n_windows)])
    else:
        pos_avg = positions
        depth_avg = depths
    
    color = '#2196F3' if seq_type == 'illumina' else '#FF5722'
    ax.fill_between(pos_avg / 1000, depth_avg, alpha=0.3, color=color)
    ax.plot(pos_avg / 1000, depth_avg, linewidth=0.5, color=color)
    
    ax.axhline(y=mean_d, color='red', linestyle='--', linewidth=0.8, alpha=0.7, label=f'Mean: {mean_d:.0f}X')
    ax.axhline(y=median_d, color='green', linestyle=':', linewidth=0.8, alpha=0.7, label=f'Median: {median_d:.0f}X')
    
    # Mark zero coverage regions
    gap_info = gap_report.get(sample, {})
    for region in gap_info.get('zero_regions', []):
        ax.axvspan(region[0]/1000, region[1]/1000, alpha=0.3, color='red', label='0X gap' if region == gap_info['zero_regions'][0] else '')
    
    ax.set_xlabel('Position (kb)', fontsize=11)
    ax.set_ylabel('Coverage Depth', fontsize=11)
    ax.set_title(f'{sample} ({seq_type.upper()}) - Mitochondrial Coverage Depth', fontsize=13, fontweight='bold')
    ax.legend(loc='upper right', fontsize=9)
    ax.set_xlim(0, positions[-1] / 1000)
    ax.set_ylim(0, None)
    
    plt.tight_layout()
    fig.savefig(os.path.join(PLOT_DIR, f'{sample}_coverage.png'), dpi=150, bbox_inches='tight')
    plt.close(fig)
    print(f"  {sample}_coverage.png")

# --- 4b. Summary plot: all samples in one figure ---
fig, axes = plt.subplots(7, 4, figsize=(24, 28), constrained_layout=True)
fig.suptitle('Mitochondrial Coverage Depth - All 28 Strains (100bp window average)', 
             fontsize=16, fontweight='bold', y=1.01)

sample_order = ['f2', 'f4', 'nn12-1', 'nn12-17'] + SRR12151_SAMPLES + SRR8699_SAMPLES

for idx, sample in enumerate(sample_order):
    row, col = idx // 4, idx % 4
    ax = axes[row][col]
    
    if sample not in all_depth_data:
        ax.text(0.5, 0.5, f'{sample}\n(no data)', ha='center', va='center', transform=ax.transAxes)
        continue
    
    positions, depths = all_depth_data[sample]
    seq_type = raw_stats.get(sample, {}).get('seq_type', 'unknown')
    mean_d = np.mean(depths)
    
    window = 200
    n_windows = len(depths) // window
    pos_avg = np.array([positions[i*window] for i in range(n_windows)])
    depth_avg = np.array([np.mean(depths[i*window:(i+1)*window]) for i in range(n_windows)])
    
    color = '#2196F3' if seq_type == 'illumina' else '#FF5722'
    ax.fill_between(pos_avg / 1000, depth_avg, alpha=0.3, color=color)
    ax.plot(pos_avg / 1000, depth_avg, linewidth=0.4, color=color)
    ax.axhline(y=mean_d, color='red', linestyle='--', linewidth=0.6, alpha=0.5)
    
    # Mark gaps
    gap_info = gap_report.get(sample, {})
    for region in gap_info.get('zero_regions', []):
        ax.axvspan(region[0]/1000, region[1]/1000, alpha=0.4, color='red')
    
    has_gap = '⚠' if gap_info.get('zero_bp', 0) > 0 else '✓'
    ax.set_title(f'{sample} ({mean_d:.0f}X) {has_gap}', fontsize=9, fontweight='bold')
    ax.set_xlim(0, positions[-1] / 1000)
    ax.set_ylim(0, None)
    ax.tick_params(labelsize=7)
    
    if row == 6:
        ax.set_xlabel('Position (kb)', fontsize=8)
    if col == 0:
        ax.set_ylabel('Depth', fontsize=8)

plt.savefig(os.path.join(PLOT_DIR, 'all_samples_coverage_summary.png'), dpi=200, bbox_inches='tight')
plt.close()
print(f"\n  all_samples_coverage_summary.png")

# ============================================================
# 5. 断裂区域详细报告
# ============================================================
print("\n" + "=" * 60)
print("断裂区域详细报告")
print("=" * 60)

gap_samples = {s: g for s, g in gap_report.items() if g['zero_bp'] > 0}
if gap_samples:
    print(f"\n有断裂(0覆盖)的菌株: {len(gap_samples)} 个\n")
    for sample, info in sorted(gap_samples.items()):
        print(f"  {sample}:")
        for start, end, length in info['zero_regions']:
            print(f"    [{start:,} - {end:,}] ({length} bp)")
        print(f"    总计: {info['zero_bp']} bp 无覆盖")
else:
    print("\n所有菌株均无断裂(0覆盖)区域 ✓")

# Low coverage report
low_samples = {s: g for s, g in gap_report.items() if g['low_bp'] > 0}
if low_samples:
    print(f"\n有低覆盖(<{LOW_COV_THRESHOLD}X)区域的菌株: {len(low_samples)} 个\n")
    for sample, info in sorted(low_samples.items()):
        total_low = info['low_bp']
        pct = total_low / info['genome_len'] * 100
        print(f"  {sample}: {total_low} bp ({pct:.3f}%) 低覆盖, {len(info['low_regions'])} 个区域")
        if len(info['low_regions']) <= 10:
            for start, end, length in info['low_regions']:
                print(f"    [{start:,} - {end:,}] ({length} bp)")

print("\n统计完成！")
print(f"结果表: {tsv_path}")
print(f"图表目录: {PLOT_DIR}/")
