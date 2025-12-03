# 从基因组中筛选CDS序列，或者是全部的基因的CDS序列
# 修改所有的起始密码子为ATG  就是那些GTG的，改为ATG
cat *.fas > all.fasta
grep -v "^>" all.fasta > xin.fasta




(base) maruichen@maruichen-X99-QD4:~/workspace/work/32.mito.mimazi$ cat mimazi.sh 
# 从基因组中筛选CDS序列，或者是全部的基因的CDS序列
# 修改所有的起始密码子为ATG  就是那些GTG的，改为ATG
cat *.fas > all.fasta
grep -v "^>" all.fasta > xin.fasta

