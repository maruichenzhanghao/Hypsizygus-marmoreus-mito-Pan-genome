#!/bin/bash

# 检查是否提供了目标目录参数
if [[ -z "$1" ]]; then
    echo "Usage: $0 <target_directory>"
    exit 1
fi

# 设置目标目录为传入的第一个参数
TARGET_DIR="$1"

# 进入目标目录
cd "$TARGET_DIR" || { echo "Error: Failed to enter directory $TARGET_DIR"; exit 1; }

# 遍历所有 .fasta 文件
for fasta_file in *.fasta; do
    # 提取文件名（不带扩展名）
    base_name=$(basename "$fasta_file" .fasta)
    
    # 定义日志文件名
    log_file="${base_name}.log"
    
    # 运行 mfannot 命令，并将输出保存到日志文件
    echo "Running mfannot on $fasta_file..."
    mfannot -g 4 --tbl --sqn "$fasta_file" > "$log_file" 2>&1
    
    # 检查命令是否成功
    if [[ $? -eq 0 ]]; then
        echo "Success: $fasta_file -> $log_file"
    else
        echo "Error: Failed to process $fasta_file"
    fi
done

echo "All files have been processed."
