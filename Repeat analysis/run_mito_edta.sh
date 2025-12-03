conda activate EDTA



EDTA.pl --genome /home/maruichen/data/02.mito/f2.fasta --species others --threads 6 --anno 1 --force 1 --evaluate 1

EDTA.pl --genome /home/maruichen/data/02.mito/f2.fasta --species others --threads 6 --anno 1 --force 1 --sensitive 1 --debug 1 --overwrite 1 --evaluate 1

cat  f2.fasta.mod.EDTA.TEanno.sum

