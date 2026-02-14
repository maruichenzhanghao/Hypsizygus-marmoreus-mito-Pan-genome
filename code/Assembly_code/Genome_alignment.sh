# 对组装结果建库
makeblastdb -in assembly.fasta -out assembly.fasta -dbtype nucl

# 比对相似序列
blastn -query ../0.cankao/fujiannonglin2020.fasta -db assembly.fasta -evalue 1e-10 -out output.txt -outfmt 6 -num_threads 20 -max_hsps 1
