#!/bin/bash

# 检查参数数量
if [ "$#" -ne 3 ]; then
    echo "Usage: $0 <input_fasta> <input_fastq> <output_dir>"
    exit 1
fi

# 获取参数
INPUT_FASTA=$1
INPUT_FASTQ=$2
OUTPUT_DIR=$3

# 创建输出目录（如果不存在）
mkdir -p "$OUTPUT_DIR"

# 进入输出目录
cd "$OUTPUT_DIR" || { echo "Failed to change directory to $OUTPUT_DIR"; exit 1; }

# 使用 minimap2 建库
minimap2 -d raw "$INPUT_FASTA"

# 使用 minimap2 进行比对
minimap2 -ax map-pb -t 4 raw "$INPUT_FASTQ" > out.asm

# 将 SAM 格式转换为 BAM 格式
samtools view -@ 6 -bS -F 4 out.asm > out.bam

# 对 BAM 文件进行排序
samtools sort -o out.sort.bam out.bam

# 对排序后的 BAM 文件进行索引
samtools index out.sort.bam

# 输出深度文档
samtools depth out.sort.bam > out.depth.txt

echo "Processing completed. Output files are in $OUTPUT_DIR"
