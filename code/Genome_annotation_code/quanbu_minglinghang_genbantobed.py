# -*- coding: utf-8 -*-

from Bio import SeqIO
import pandas as pd
import glob
import os
import sys

def extract_genes_from_genbank(gb_file):
    data = []
    for record in SeqIO.parse(gb_file, "genbank"):
        for feature in record.features:
            if feature.type == "gene":
                gene_name = feature.qualifiers.get("gene", [""])[0]
                if not gene_name:
                    gene_name = feature.qualifiers.get("locus_tag", [""])[0]
                data.append({
                    "molecule": record.id,
                    "gene": gene_name,
                    "start": int(feature.location.start) + 1,
                    "end": int(feature.location.end),
                    "strand": "forward" if feature.location.strand > 0 else "reverse",
                    "orientation": 1 if feature.location.strand < 0 else 0
                })
    return pd.DataFrame(data)

def process_genbank_folder(input_dir=".", output_dir="output"):
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    
    genbank_files = glob.glob(os.path.join(input_dir, "*.gb")) + glob.glob(os.path.join(input_dir, "*.gbk"))
    
    for gb_file in genbank_files:
        df = extract_genes_from_genbank(gb_file)
        base_name = os.path.basename(gb_file).rsplit(".", 1)[0]
        output_csv = os.path.join(output_dir, f"{base_name}_genes.csv")
        df.to_csv(output_csv, index=False)
        print(f"Processed: {gb_file} -> {output_csv}")

if __name__ == "__main__":
    input_dir = sys.argv[1] if len(sys.argv) > 1 else "."
    output_dir = sys.argv[2] if len(sys.argv) > 2 else "output"
    process_genbank_folder(input_dir, output_dir)
    print("Done!")  # 替换原来的中文输出
