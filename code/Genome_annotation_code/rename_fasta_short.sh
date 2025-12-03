#!/bin/bash

# 遍历当前目录下的所有 .fasta 文件
for file in *.fasta; do
    # 获取文件名（去掉 .fasta 后缀），并提取 .K105 之前的部分
    name=$(basename "$file" .fasta | cut -d'.' -f1)
    
    # 使用 sed 替换第一行
    sed -i "1s/^>.*$/>$name/" "$file"
done
