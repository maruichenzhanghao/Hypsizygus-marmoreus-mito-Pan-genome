#!/bin/bash

# 遍历当前目录下所有的 .sqn 文件
for sqn_file in *.sqn; do
    # 提取文件名前缀（去掉 .fasta.new.sqn 部分）
    prefix=$(basename "$sqn_file" .fasta.new.sqn)
    
    # 定义输出文件名
    output_file="${prefix}.gb"
    
    # 使用 asn2gb 命令进行转换
    asn2gb -i "$sqn_file" -o "$output_file" -f b -h 0 -l logfile.txt
    
    # 打印信息以便跟踪进度
    echo "Converted $sqn_file to $output_file"
done
