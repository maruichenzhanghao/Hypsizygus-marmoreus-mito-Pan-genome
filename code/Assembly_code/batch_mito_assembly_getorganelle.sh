#!/bin/bash

# 检查是否提供了输出目录作为参数
if [ -z "$1" ]; then
    echo "Usage: $0 <output_directory>"
    exit 1
fi

# 获取输出目录路径
OUTPUT_DIR="$1"

# 创建输出目录（如果不存在）
mkdir -p "$OUTPUT_DIR"

# 定义参考序列和其他固定参数
REFERENCE="all.fasta"
KMER_SIZES="21,45,65,85,105"
FUNGUS_MT_FLAG="-F fungus_mt"
READ_REPETITIONS="-R 10"
CONTINUE_OPTION="--continue"  # 可选：使用 '--overwrite' 如果需要覆盖
TIMEOUT_DURATION="30m"       # 设置超时时间为30分钟

# 获取当前目录下的所有 _1.fastq.gz 文件列表
FASTQ_FILES_1=($(ls *_1.fastq.gz 2>/dev/null))

# 检查是否有 _1.fastq.gz 文件存在
if [ ${#FASTQ_FILES_1[@]} -eq 0 ]; then
    echo "No _1.fastq.gz files found in the current directory."
    exit 1
fi

# 遍历每个 _1.fastq.gz 文件并找到对应的 _2.fastq.gz 文件
for FILE_1 in "${FASTQ_FILES_1[@]}"; do
    # 获取基础文件名（不包括 _1.fastq.gz）
    BASENAME="${FILE_1%%_1.fastq.gz}"
    
    # 构建对应的 _2.fastq.gz 文件名
    FILE_2="${BASENAME}_2.fastq.gz"
    
    # 检查对应的 _2.fastq.gz 文件是否存在
    if [ ! -f "$FILE_2" ]; then
        echo "Warning: Corresponding file for $FILE_1 not found: $FILE_2. Skipping..."
        continue
    fi
    
    # 构建输出子目录名
    OUTPUT_SUBDIR="$OUTPUT_DIR/${BASENAME}_fungus_mt_out"
    
    # 创建输出子目录（如果不存在）
    mkdir -p "$OUTPUT_SUBDIR"
    
    # 构建并执行 get_organelle_from_reads.py 命令
    echo "Starting processing for $BASENAME..."

    # 使用 timeout 命令限制运行时间
    timeout $TIMEOUT_DURATION get_organelle_from_reads.py \
        -1 "$FILE_1" \
        -2 "$FILE_2" \
        $READ_REPETITIONS \
        -s "$REFERENCE" \
        -k "$KMER_SIZES" \
        $FUNGUS_MT_FLAG \
        $CONTINUE_OPTION \
        -o "$OUTPUT_SUBDIR"
    
    # 检查命令是否成功执行
    if [ $? -eq 0 ]; then
        echo "Successfully processed $BASENAME."
    elif [ $? -eq 124 ]; then
        echo "Error: Processing of $BASENAME exceeded the time limit and was terminated."
    else
        echo "Error: Failed to process $BASENAME."
        read -p "Press Enter to continue with the next sample or Ctrl+C to abort..."
    fi
    
done

echo "All samples have been processed."
