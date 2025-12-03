seqkit sample -p 0.2 -s 11 f2.fastq.gz > f2.sample.fq
flye --pacbio-hifi /home/maruichen/data/2.3dai/F2/f2.sample.fq --out-dir /home/maruichen/workspace/1zz.mito/1.3daizz/1.flye --threads 16
