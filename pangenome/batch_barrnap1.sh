#!/bin/bash

##############################################################################
# 功能：批量对指定目录下的.fasta文件运行barrnap命令，生成rRNA序列和GFF3注释
# 使用方式：./batch_barrnap.sh [输入目录] [输出目录（可选）]
# 示例1：./batch_barrnap.sh ~/input_files ~/output_results  # 指定输出目录
# 示例2：./batch_barrnap.sh ~/input_files                    # 输出到当前目录
##############################################################################

# 1. 检查输入参数（至少需要输入目录）
if [ $# -lt 1 ]; then
    echo "ERROR: 请传入输入目录路径！"
    echo "使用示例："
    echo "  指定输出目录：./batch_barrnap.sh 输入目录 输出目录"
    echo "  默认输出目录：./batch_barrnap.sh 输入目录"
    exit 1
fi

# 2. 定义输入目录并检查是否存在
INPUT_DIR="$1"
if [ ! -d "$INPUT_DIR" ]; then
    echo "ERROR: 输入目录不存在！目录路径：$INPUT_DIR"
    exit 1
fi

# 3. 定义输出目录（若未指定则使用当前目录）
if [ $# -eq 2 ]; then
    OUTPUT_DIR="$2"
else
    OUTPUT_DIR=$(pwd)  # 当前目录作为默认输出目录
fi

# 4. 创建输出目录（若不存在）
mkdir -p "$OUTPUT_DIR" || { echo "ERROR: 无法创建输出目录 $OUTPUT_DIR"; exit 1; }

# 5. 查找输入目录下所有.fasta文件（仅当前目录，不递归子目录）
FASTA_FILES=$(find "$INPUT_DIR" -maxdepth 1 -type f -name "*.fasta")

# 6. 检查是否有.fasta文件
if [ -z "$FASTA_FILES" ]; then
    echo "WARNING: 输入目录下未找到.fasta文件！目录路径：$INPUT_DIR"
    exit 0
fi

# 7. 批量运行barrnap命令
echo "===== 开始批量处理 ====="
echo "输入目录：$INPUT_DIR"
echo "输出目录：$OUTPUT_DIR"
echo "======================"

for file in $FASTA_FILES; do
    # 获取文件名（不含路径和扩展名）
    filename=$(basename "$file" .fasta)
    # 定义输出文件路径
    outseq="$OUTPUT_DIR/${filename}_rRNA.fasta"
    outgff="$OUTPUT_DIR/${filename}_rRNA.gff3"

    echo -e "\n----- 正在处理文件：$filename.fasta -----"
    echo "输出rRNA序列：$outseq"
    echo "输出GFF3注释：$outgff"

    # 运行barrnap（已移除--incseq选项）
    barrnap --kingdom euk \
            --threads 10 \
            --outseq "$outseq" \
            "$file" > "$outgff"

    # 检查命令执行结果
    if [ $? -eq 0 ]; then
        echo "处理完成：$filename.fasta"
    else
        echo "ERROR: 处理失败！文件：$filename.fasta"
    fi
done

echo -e "\n===== 批量处理结束！所有输出文件已保存至：$OUTPUT_DIR ====="

