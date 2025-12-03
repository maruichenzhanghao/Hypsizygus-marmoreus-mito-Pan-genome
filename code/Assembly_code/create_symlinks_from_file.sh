#!/bin/bash

# 检查是否提供了输入文件作为参数
if [ -z "$1" ]; then
    echo "Usage: $0 <path_to_file_with_directories>"
    exit 1
fi

# 获取输入文件路径
INPUT_FILE="$1"

# 检查文件是否存在且可读
if [ ! -f "$INPUT_FILE" ] || [ ! -r "$INPUT_FILE" ]; then
    echo "Error: File '$INPUT_FILE' does not exist or is not readable."
    exit 1
fi

# 遍历输入文件中的每一行（每个目录）
while IFS= read -r DIR; do
    # 忽略空行和以#开头的注释行
    if [[ -n "$DIR" && ! "$DIR" =~ ^# ]]; then
        # 查找该目录下的所有 .fastq.gz 文件
        find "$DIR" -type f -name "*.fastq.gz" | while read -r FILE; do
            # 获取文件名（不包括路径）
            BASENAME=$(basename "$FILE")
            
            # 检查是否已经在当前目录存在同名文件或链接
            if [ ! -e "./$BASENAME" ]; then
                # 创建软链接
                ln -s "$FILE" .
                echo "Created symlink for $FILE"
            else
                echo "File or symlink './$BASENAME' already exists, skipping."
            fi
        done
    fi
done < "$INPUT_FILE"
