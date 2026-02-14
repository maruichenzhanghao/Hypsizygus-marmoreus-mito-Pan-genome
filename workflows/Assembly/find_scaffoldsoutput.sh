#!/bin/bash

# 定义目标输出目录
OUTPUT_DIR="scaffoldsoutput"

# 创建输出目录，如果它不存在的话
mkdir -p "$OUTPUT_DIR"

# 遍历所有_fungus_mt_out结尾的目录
for dir in *_fungus_mt_out; do
    # 检查是否为目录
    if [ -d "$dir" ]; then
        # 获取目录名作为前缀
        prefix="${dir%_fungus_mt_out}"

        # 查找并处理符合模式的文件
        find "$dir" -type f -name '*fungus_mt.K*.scaffolds.graph*.path_sequence.fasta' | while read -r file; do
            # 获取文件的基本名和目录
            base_name=$(basename "$file")
            dir_name=$(dirname "$file")

            # 构造新的文件名
            new_name="${prefix}${base_name}"

            # 复制并重命名文件到输出目录
            cp "$file" "$OUTPUT_DIR/$new_name"
            
            echo "Copied and renamed: $file -> $OUTPUT_DIR/$new_name"
        done
    fi
done

echo "操作完成！所有匹配的文件已复制并重命名至 $OUTPUT_DIR 目录。"

