import os
from Bio import SeqIO
import uuid

# 输入和输出目录
input_dir = os.path.expanduser("~/workspace/01.mito/02.gb/1.zuizhong.8.19")  # 使用os.path.expanduser处理~
output_dir = os.path.expanduser("~/workspace/01.mito/02.gb/fasta_output")  # 输出FASTA目录

# 确保输出目录存在
os.makedirs(output_dir, exist_ok=True)

# 获取所有GB文件
gb_files = [f for f in os.listdir(input_dir) if f.endswith('.gb')]

for gb_file in gb_files:
    input_path = os.path.join(input_dir, gb_file)
    output_fasta = os.path.join(output_dir, f"{gb_file.replace('.gb', '')}.fasta")
    
    # 初始化FASTA内容
    fasta_records = []
    
    # 解析GenBank文件
    for record in SeqIO.parse(input_path, "genbank"):
        for feature in record.features:
            if feature.type == "CDS":
                # 获取蛋白序列
                protein_seq = feature.qualifiers.get('translation', [None])[0]
                if not protein_seq:
                    continue  # 跳过没有翻译序列的CDS
                
                # 获取基因名或其他标识
                gene_name = feature.qualifiers.get('gene', ['unknown'])[0]
                product = feature.qualifiers.get('product', ['unknown'])[0]
                # 生成唯一ID（避免重复）
                seq_id = f"{gb_file.replace('.gb', '')}_{gene_name}_{product}_{uuid.uuid4().hex[:8]}"
                
                # 创建FASTA记录
                fasta_record = f">{seq_id}\n{protein_seq}\n"
                fasta_records.append(fasta_record)
    
    # 写入FASTA文件
    with open(output_fasta, 'w') as f:
        f.writelines(fasta_records)
    
    print(f"Generated FASTA file: {output_fasta}")

print("All GB files processed. Ready for OrthoFinder analysis.")

# 下一步：运行OrthoFinder
print("""
To perform pangenome analysis with OrthoFinder:
1. Ensure OrthoFinder is installed (`pip install orthofinder` or from GitHub).
2. Move all generated FASTA files to a single directory (e.g., fasta_output).
3. Run OrthoFinder:
   orthofinder -f ~/workspace/01.mito/02.gb/fasta_output
4. OrthoFinder will output orthogroups in the results directory.

To visualize core/shell/cloud genes:
1. Use the Orthogroups.GeneCount.tsv file from OrthoFinder output.
2. Parse it with Python/R to classify genes:
   - Core: present in all 31 strains (or >95% threshold).
   - Shell: present in 2-30 strains.
   - Cloud: present in 1 strain.
3. Visualize with R (ggplot2/UpSetR) or Python (matplotlib/seaborn).
""")
