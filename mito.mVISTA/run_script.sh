#!/bin/bash

# 遍历当前目录下所有的.gb文件
for gbfile in *.gb; do
    # 执行你的Python脚本
    python get_mVISTA_annotation_file_from_genbank_1.py -i "$gbfile"
done
