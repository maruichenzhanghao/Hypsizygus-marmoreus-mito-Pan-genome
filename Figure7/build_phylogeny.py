#!/usr/bin/env python3
"""
白玉菇线粒体基因组系统发育分析流程
=====================================
从31个GenBank文件中提取15个核心蛋白编码基因（排除dpo），
逐基因比对后串联(concatenation)，使用IQ-TREE构建系统发育树。

15个核心基因：
cox1, cox2, cox3, cob, nad1, nad2, nad3, nad4, nad4L, nad5, nad6, 
atp6, atp8, atp9, rps3

作者: 自动生成
日期: 2026-02-11
"""

import os
import sys
import subprocess
from pathlib import Path
from collections import defaultdict

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

# ==================== 配置 ====================
GB_DIR = Path("/home/maxinxin/workspace/01.mito/02.gb/1.gb/1.zuizhong.8.19")
OUT_DIR = Path("/home/maxinxin/workspace/01.mito/03.pangenome.cactus/phylogeny")
OUT_DIR.mkdir(parents=True, exist_ok=True)

# 15个核心蛋白编码基因（排除dpo）
CORE_GENES = [
    "cox1", "cox2", "cox3", "cob",
    "nad1", "nad2", "nad3", "nad4", "nad4L", "nad5", "nad6",
    "atp6", "atp8", "atp9", "rps3"
]

# 31个样本
SAMPLES = [
    "f2", "f4", "MF133443.1", "MH382825.1",
    "nn12-1", "nn12-17",
    "SRR12151860", "SRR12151871", "SRR12151875", "SRR12151883",
    "SRR7874787",
    "SRR8699796", "SRR8699797", "SRR8699800", "SRR8699801",
    "SRR8699802", "SRR8699803", "SRR8699804", "SRR8699805",
    "SRR8699808", "SRR8699809", "SRR8699811", "SRR8699813",
    "SRR8699814", "SRR8699815", "SRR8699816", "SRR8699817",
    "SRR8699833", "SRR8699834", "SRR8699835", "SRR8699837"
]


def extract_cds_sequence(gb_file, gene_name):
    """
    从GenBank文件中提取指定基因的CDS核苷酸序列。
    处理含内含子的基因（join()位置），只提取外显子拼接后的编码序列。
    """
    for record in SeqIO.parse(gb_file, "genbank"):
        for feature in record.features:
            if feature.type == "CDS":
                gene = feature.qualifiers.get("gene", [""])[0]
                if gene.lower() == gene_name.lower():
                    # 提取拼接后的核苷酸序列
                    seq = feature.location.extract(record.seq)
                    return str(seq)
    return None


def step1_extract_genes():
    """步骤1: 从GenBank文件中提取15个核心基因的CDS序列"""
    print("=" * 60)
    print("步骤1: 提取15个核心基因的CDS序列")
    print("=" * 60)
    
    gene_dir = OUT_DIR / "01_gene_sequences"
    gene_dir.mkdir(exist_ok=True)
    
    # 统计矩阵
    gene_matrix = {}  # {gene: {sample: seq}}
    missing = []
    
    for gene in CORE_GENES:
        gene_matrix[gene] = {}
    
    for sample in SAMPLES:
        gb_file = GB_DIR / f"{sample}.gb"
        if not gb_file.exists():
            print(f"  [警告] 找不到文件: {gb_file}")
            continue
        
        for gene in CORE_GENES:
            seq = extract_cds_sequence(gb_file, gene)
            if seq:
                gene_matrix[gene][sample] = seq
            else:
                missing.append((sample, gene))
                print(f"  [缺失] {sample} 缺少 {gene}")
    
    # 写出每个基因的FASTA文件
    for gene in CORE_GENES:
        fasta_file = gene_dir / f"{gene}.fasta"
        with open(fasta_file, "w") as f:
            for sample in SAMPLES:
                if sample in gene_matrix[gene]:
                    f.write(f">{sample}\n{gene_matrix[gene][sample]}\n")
        
        n_seqs = len(gene_matrix[gene])
        print(f"  {gene}: {n_seqs}/31 个样本提取成功, 序列写入 {fasta_file.name}")
    
    # 打印统计
    print(f"\n  总计缺失: {len(missing)} 个基因-样本组合")
    if missing:
        for s, g in missing:
            print(f"    - {s}: {g}")
    
    return gene_matrix


def step2_align_genes():
    """步骤2: 使用MAFFT对每个基因进行多序列比对"""
    print("\n" + "=" * 60)
    print("步骤2: MAFFT多序列比对（逐基因）")
    print("=" * 60)
    
    gene_dir = OUT_DIR / "01_gene_sequences"
    aln_dir = OUT_DIR / "02_alignments"
    aln_dir.mkdir(exist_ok=True)
    
    for gene in CORE_GENES:
        input_fasta = gene_dir / f"{gene}.fasta"
        output_aln = aln_dir / f"{gene}_aligned.fasta"
        
        if not input_fasta.exists():
            print(f"  [跳过] {gene}: 输入文件不存在")
            continue
        
        # 使用MAFFT进行比对（--auto自动选择最佳策略）
        cmd = [
            "mafft", "--auto", "--adjustdirectionaccurately",
            "--thread", "4",
            str(input_fasta)
        ]
        
        try:
            result = subprocess.run(
                cmd, capture_output=True, text=True, timeout=300
            )
            
            # MAFFT输出到stdout
            with open(output_aln, "w") as f:
                f.write(result.stdout)
            
            # 检查结果
            aligned_seqs = list(SeqIO.parse(output_aln, "fasta"))
            if aligned_seqs:
                aln_len = len(aligned_seqs[0].seq)
                print(f"  {gene}: 比对完成, {len(aligned_seqs)}条序列, 比对长度={aln_len} bp")
            else:
                print(f"  [错误] {gene}: 比对结果为空")
        except subprocess.TimeoutExpired:
            print(f"  [超时] {gene}: MAFFT比对超时")
        except Exception as e:
            print(f"  [错误] {gene}: {e}")


def step3_concatenate():
    """步骤3: 串联所有基因比对结果，生成partition文件"""
    print("\n" + "=" * 60)
    print("步骤3: 串联比对序列 + 生成partition文件")
    print("=" * 60)
    
    aln_dir = OUT_DIR / "02_alignments"
    concat_dir = OUT_DIR / "03_concatenated"
    concat_dir.mkdir(exist_ok=True)
    
    # 读取所有比对
    all_alignments = {}  # {gene: {sample: aligned_seq}}
    gene_lengths = {}
    
    for gene in CORE_GENES:
        aln_file = aln_dir / f"{gene}_aligned.fasta"
        if not aln_file.exists():
            print(f"  [跳过] {gene}: 比对文件不存在")
            continue
        
        seqs = {}
        for record in SeqIO.parse(aln_file, "fasta"):
            # MAFFT --adjustdirectionaccurately可能在ID前加_R_前缀表示反向
            sample_id = record.id.replace("_R_", "")
            seqs[sample_id] = str(record.seq)
        
        if seqs:
            aln_len = len(list(seqs.values())[0])
            all_alignments[gene] = seqs
            gene_lengths[gene] = aln_len
            print(f"  {gene}: 读取{len(seqs)}条序列, 长度={aln_len}")
    
    # 确定所有样本集合
    all_samples = set()
    for gene_seqs in all_alignments.values():
        all_samples.update(gene_seqs.keys())
    all_samples = sorted(all_samples)
    
    print(f"\n  总样本数: {len(all_samples)}")
    print(f"  总基因数: {len(all_alignments)}")
    
    # 串联
    concatenated = {}
    for sample in all_samples:
        concat_seq = ""
        for gene in CORE_GENES:
            if gene in all_alignments:
                if sample in all_alignments[gene]:
                    concat_seq += all_alignments[gene][sample]
                else:
                    # 用gap填充缺失基因
                    concat_seq += "-" * gene_lengths[gene]
                    print(f"  [填充] {sample} 缺少 {gene}，用gap填充")
            
        concatenated[sample] = concat_seq
    
    total_len = sum(gene_lengths.values())
    print(f"\n  串联总长度: {total_len} bp")
    
    # 写出串联FASTA
    concat_fasta = concat_dir / "15genes_concatenated.fasta"
    with open(concat_fasta, "w") as f:
        for sample in all_samples:
            f.write(f">{sample}\n")
            # 每行80个字符
            seq = concatenated[sample]
            for i in range(0, len(seq), 80):
                f.write(seq[i:i+80] + "\n")
    print(f"  串联序列写入: {concat_fasta}")
    
    # 生成partition文件（RAxML/IQ-TREE格式）
    partition_file = concat_dir / "partition.txt"
    pos = 1
    with open(partition_file, "w") as f:
        for gene in CORE_GENES:
            if gene in gene_lengths:
                end = pos + gene_lengths[gene] - 1
                f.write(f"DNA, {gene} = {pos}-{end}\n")
                pos = end + 1
    print(f"  Partition文件写入: {partition_file}")
    
    # 生成NEXUS格式partition（IQ-TREE也支持）
    nexus_part = concat_dir / "partition.nex"
    with open(nexus_part, "w") as f:
        f.write("#nexus\nbegin sets;\n")
        pos = 1
        for gene in CORE_GENES:
            if gene in gene_lengths:
                end = pos + gene_lengths[gene] - 1
                f.write(f"  charset {gene} = {pos}-{end};\n")
                pos = end + 1
        # charpartition
        gene_list = [g for g in CORE_GENES if g in gene_lengths]
        f.write("  charpartition mine = " + 
                ", ".join(f"GTR+G:{g}" for g in gene_list) + ";\n")
        f.write("end;\n")
    print(f"  NEXUS partition写入: {nexus_part}")
    
    return concat_fasta, partition_file


def step4_build_tree(concat_fasta, partition_file):
    """步骤4: 使用IQ-TREE构建系统发育树"""
    print("\n" + "=" * 60)
    print("步骤4: IQ-TREE系统发育分析")
    print("=" * 60)
    
    tree_dir = OUT_DIR / "04_phylogeny"
    tree_dir.mkdir(exist_ok=True)
    
    prefix = tree_dir / "mito_15genes"
    
    # ---- 4a: 串联无分区 (简单模型) ----
    print("\n  [4a] 运行IQ-TREE（串联, ModelFinder自动选模型）...")
    cmd_simple = [
        "iqtree",
        "-s", str(concat_fasta),
        "--prefix", str(prefix) + "_concat",
        "-m", "MFP",            # ModelFinder Plus: 自动选最佳模型
        "-B", "1000",            # UFBoot2 1000次
        "--alrt", "1000",        # SH-aLRT 1000次
        "-T", "AUTO",            # 自动线程
        "--redo"                 # 覆盖已有结果
    ]
    
    print(f"  命令: {' '.join(cmd_simple)}")
    result = subprocess.run(cmd_simple, capture_output=True, text=True, timeout=3600)
    if result.returncode == 0:
        print(f"  ✓ 串联分析完成")
    else:
        print(f"  ✗ 错误: {result.stderr[-500:]}")
    
    # ---- 4b: 分区模型 ----
    print("\n  [4b] 运行IQ-TREE（分区模型, 每个基因独立选模型）...")
    cmd_partition = [
        "iqtree",
        "-s", str(concat_fasta),
        "-p", str(partition_file),  # 分区文件
        "--prefix", str(prefix) + "_partition",
        "-m", "MFP",
        "-B", "1000",
        "--alrt", "1000",
        "-T", "AUTO",
        "--redo"
    ]
    
    print(f"  命令: {' '.join(cmd_partition)}")
    result = subprocess.run(cmd_partition, capture_output=True, text=True, timeout=3600)
    if result.returncode == 0:
        print(f"  ✓ 分区分析完成")
    else:
        print(f"  ✗ 错误: {result.stderr[-500:]}")
    
    # 打印结果文件
    print("\n  结果文件:")
    for ext in [".treefile", ".iqtree", ".log", ".contree"]:
        for tag in ["_concat", "_partition"]:
            f = Path(str(prefix) + tag + ext)
            if f.exists():
                print(f"    {f.name} ({f.stat().st_size:,} bytes)")


def step5_summary():
    """步骤5: 输出最终树文件路径"""
    print("\n" + "=" * 60)
    print("步骤5: 分析完成总结")
    print("=" * 60)
    
    tree_dir = OUT_DIR / "04_phylogeny"
    
    print("\n  关键结果文件:")
    print(f"  ├── 基因序列:    {OUT_DIR / '01_gene_sequences/'}")
    print(f"  ├── 比对文件:    {OUT_DIR / '02_alignments/'}")
    print(f"  ├── 串联序列:    {OUT_DIR / '03_concatenated/'}")
    print(f"  └── 系统发育树:  {tree_dir}/")
    
    # 读取并打印最佳树
    for tag in ["_concat", "_partition"]:
        tree_file = tree_dir / f"mito_15genes{tag}.treefile"
        if tree_file.exists():
            with open(tree_file) as f:
                tree = f.read().strip()
            print(f"\n  {'串联' if 'concat' in tag else '分区'}模型最优树 (Newick):")
            print(f"  {tree[:200]}...")
            
            # 读取IQ-TREE报告中的模型信息
            iqtree_file = tree_file.with_suffix(".iqtree")
            if iqtree_file.exists():
                with open(iqtree_file) as f:
                    for line in f:
                        if "Best-fit model" in line:
                            print(f"  最优模型: {line.strip()}")
                            break
    
    print("\n  提示:")
    print("  - .treefile 文件可直接导入 iTOL (https://itol.embl.de/) 或 FigTree 可视化")
    print("  - .contree 文件为一致树(consensus tree)")
    print("  - .iqtree 文件含详细的模型选择和bootstrap统计")
    print("  - 分区模型通常优于串联模型，建议优先使用分区结果")


def main():
    print("白玉菇线粒体基因组 15核心基因 系统发育分析")
    print(f"GenBank目录: {GB_DIR}")
    print(f"输出目录: {OUT_DIR}")
    print(f"核心基因 ({len(CORE_GENES)}个): {', '.join(CORE_GENES)}")
    print(f"样本数: {len(SAMPLES)}")
    print()
    
    # Step 1: 提取基因
    gene_matrix = step1_extract_genes()
    
    # Step 2: MAFFT比对
    step2_align_genes()
    
    # Step 3: 串联
    concat_fasta, partition_file = step3_concatenate()
    
    # Step 4: IQ-TREE建树
    step4_build_tree(concat_fasta, partition_file)
    
    # Step 5: 总结
    step5_summary()


if __name__ == "__main__":
    main()
