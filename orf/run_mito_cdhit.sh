#!/bin/bash

# $1 是传入的第一个参数，代表输入文件的路径
input_file="$1"  
output_file="clustered_orfs"

# 确保输入文件路径不是空的
if [ -z "$input_file" ]; then
  echo "错误：没有提供输入文件路径"
  exit 1
fi

# 运行cd-hit命令
/home/maruichen/software/6.annotation/08.cdhit/cd-hit-v4.8.1-2019-0228/cd-hit -i "$input_file" -o "$output_file" -c 0.8 -s 0.8 -n 5 -M 16000 -T 16
