#!/bin/bash

# 检查是否提供了正确的参数数量
if [ "$#" -ne 2 ]; then
    echo "Usage: $0 <search_directory> <output_file>"
    exit 1
fi

# 获取命令行参数
SEARCH_DIR="$1"
OUTPUT_FILE="$2"

# 检查提供的目录是否存在且可读
if [ ! -d "$SEARCH_DIR" ]; then
    echo "Error: Directory '$SEARCH_DIR' does not exist or is not accessible."
    exit 1
fi

# 清空或创建输出文件
> "$OUTPUT_FILE"

# 查找并处理文件
find "$SEARCH_DIR" -type f -name '*complete*' | while read -r file; do
    # 获取文件所在的目录
    dir=$(dirname "$file")
    # 将目录添加到输出文件中（如果该目录还没有被记录）
    if ! grep -Fxq "$dir" "$OUTPUT_FILE"; then
        echo "$dir" >> "$OUTPUT_FILE"
    fi
done

echo "完成！结果已保存至 $OUTPUT_FILE"
