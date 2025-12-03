#!/bin/bash

# 定义工作目录
WORK_DIR="/home/maruichen/workspace/work/13.mito.orf"

# 进入工作目录
cd "$WORK_DIR" || exit

# 遍历目录下所有的 .fas 文件
for file in *.fas; do
  # 创建临时文件
  temp_file="${file}.tmp"

  # 处理文件内容
  while IFS= read -r line; do
    if [[ $line == \>* ]]; then
      # 如果是序列标识行，则修改标识符，添加 '_orfXXX' 后缀
      id="${line#*>}"  # 去掉前面的 '>'
      echo ">${id}_orf${file%.fas}" >> "$temp_file"
    else
      # 否则直接写入序列数据
      echo "$line" >> "$temp_file"
    fi
  done < "$file"

  # 替换原文件为修改后的文件
  mv "$temp_file" "$file"
done

echo "All files have been processed."
