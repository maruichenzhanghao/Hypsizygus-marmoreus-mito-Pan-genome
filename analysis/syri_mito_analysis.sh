#!/bin/bash
set -euo pipefail  # 遇到错误立即退出，增强脚本稳定性

# -------------------------- 配置参数（可根据需要修改） --------------------------
REF="f2.fasta"  # 参考基因组文件名（确保在当前目录）
MUMMER_SIF="/home/maxinxin/software/01.singlity/mummer-alpine.sif"  # singularity镜像路径
BCFTOOLS_PATH="/home/maxinxin/software/bcftools-1.22"  # bcftools安装路径
MUMMER_PATH="/home/maxinxin/software/mummer4.01/bin"  # MUMmer4路径
CONDA_ENV="syri_env"  # syri所在的conda环境名
# ------------------------------------------------------------------------------

# 检查参考基因组是否存在
if [ ! -f "$REF" ]; then
    echo "错误：参考基因组文件 $REF 不存在！"
    exit 1
fi

# 提取参考基因组名称（用于输出文件前缀）
REF_NAME=$(basename "$REF" .fasta)

# 获取所有待分析的查询基因组（排除参考基因组和非FASTA文件）
QUERY_FILES=$(ls *.fasta | grep -v "^${REF_NAME}\.fasta$" | sort)

# 检查是否有查询文件
if [ -z "$QUERY_FILES" ]; then
    echo "警告：未找到待分析的查询基因组（仅排除了参考基因组 $REF）"
    exit 1
fi

# 激活conda环境（确保syri可用）
echo "激活conda环境 $CONDA_ENV..."
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate "$CONDA_ENV" || {
    echo "错误：激活conda环境 $CONDA_ENV 失败！"
    exit 1
}

# 添加MUMmer路径到环境变量
export PATH="$MUMMER_PATH:$PATH"
echo "MUMmer4 路径已添加：$MUMMER_PATH"

# 循环处理每个查询基因组
for QUERY in $QUERY_FILES; do
    # 提取查询基因组名称（去除.fasta后缀）
    QUERY_NAME=$(basename "$QUERY" .fasta)
    echo -e "\n======================================"
    echo "开始处理样本：$QUERY_NAME"
    echo "参考基因组：$REF_NAME"
    echo "======================================"

    # 创建输出目录（避免文件冲突）
    OUT_DIR="${REF_NAME}_vs_${QUERY_NAME}_syri_output"
    mkdir -p "$OUT_DIR"
    echo "输出目录：$OUT_DIR"

    # 1. 使用nucmer进行序列比对（ singularity执行mummer工具）
    echo "步骤1/6：运行nucmer比对..."
    singularity exec "$MUMMER_SIF" nucmer \
        --maxmatch \
        -l 20 \
        -b 200 \
        -c 65 \
        -p "$OUT_DIR/${REF_NAME}_vs_${QUERY_NAME}" \
        "$REF" \
        "$QUERY"

    # 2. 过滤delta文件
    echo "步骤2/6：过滤delta文件..."
    singularity exec "$MUMMER_SIF" delta-filter \
        -r -q \
        "$OUT_DIR/${REF_NAME}_vs_${QUERY_NAME}.delta" \
        > "$OUT_DIR/${REF_NAME}_vs_${QUERY_NAME}.filtered.delta"

    # 3. 生成坐标文件
    echo "步骤3/6：生成坐标文件..."
    singularity exec "$MUMMER_SIF" show-coords \
        -THrd \
        "$OUT_DIR/${REF_NAME}_vs_${QUERY_NAME}.filtered.delta" \
        > "$OUT_DIR/${REF_NAME}_vs_${QUERY_NAME}.filtered.coords"

    # 4. 运行syri检测结构变异
    echo "步骤4/6：运行syri分析..."
    syri \
        -c "$OUT_DIR/${REF_NAME}_vs_${QUERY_NAME}.filtered.coords" \
        -d "$OUT_DIR/${REF_NAME}_vs_${QUERY_NAME}.filtered.delta" \
        -r "$REF" \
        -q "$QUERY" \
        -F T \
        --dir "$OUT_DIR"

    # 5. 使用bcftools统计VCF并绘图
    echo "步骤5/6：生成变异统计结果..."
    export PATH="$BCFTOOLS_PATH:$PATH"  # 添加bcftools路径
    bcftools stats "$OUT_DIR/syri.vcf" > "$OUT_DIR/syri.stats.txt"
    "$BCFTOOLS_PATH/misc/plot-vcfstats" -p "$OUT_DIR/stats_plots" "$OUT_DIR/syri.stats.txt"

    # 6. 过滤VCF，保留插入（INS）和缺失（DEL）变异
    echo "步骤6/6：过滤VCF文件..."
    bcftools view -i 'ALT ~ "<INS>" || ALT ~ "<DEL>"' "$OUT_DIR/syri.vcf" -o "$OUT_DIR/syri.filtered.vcf"

    echo -e "\n样本 $QUERY_NAME 处理完成！结果保存至：$OUT_DIR"
done

echo -e "\n======================================"
echo "所有样本处理完毕！"
echo "参考基因组：$REF_NAME"
echo "分析结果均保存在各样本对应的syri_output目录中"
echo "======================================"

# 退出conda环境
conda deactivate
