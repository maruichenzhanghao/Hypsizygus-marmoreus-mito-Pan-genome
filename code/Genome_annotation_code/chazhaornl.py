from Bio import SeqIO

def find_rrna_genes(gbff_file):
    # 解析 GenBank 文件
    record = SeqIO.read(gbff_file, "genbank")
    
    # 遍历所有基因特征
    for feature in record.features:
        if feature.type == "gene":
            # 检查基因名称是否包含 rRNA 相关关键词
            gene_name = None
            for key in feature.qualifiers.keys():
                if key == "gene":
                    gene_name = feature.qualifiers["gene"][0].lower()  # 提取基因名并转为小写
                    break
            
            # 匹配 rnl、rrnL、16S 等关键词
            if gene_name and ("rnl" in gene_name or "rrnl" in gene_name or "16s" in gene_name):
                print(f"Found rRNA gene: {gene_name}")
                print(f"Location: {feature.location}")
                print(f"Qualifiers: {feature.qualifiers}\n")

# 使用示例
find_rrna_genes("f2-Reverse.gbff")
