#!/bin/bash

# 输出文件
output_file="mito_samples.txt"

# 清空输出文件（如果已存在）
> "$output_file"

# 遍历当前目录中的所有 .fasta 文件
for file in *.fasta; do
    # 获取文件名（去掉 .fasta 后缀）
    sample_name="${file%.fasta}"
    # 获取文件的绝对路径
    file_path="$(pwd)/$file"
    # 将样本名称和文件路径写入输出文件
    echo -e "$sample_name\t$file_path" >> "$output_file"
done

echo "输入文件已生成：$output_file"
