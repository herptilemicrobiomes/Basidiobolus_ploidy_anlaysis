#!/usr/bin/bash -l
#SBATCH -p short -c 8 --mem 16gb  --out logs/index.log

CPU=4
module load samtools
module load bwa-mem2

module load parallel

parallel -j $CPU samtools faidx {} ::: $(ls genome/*.fasta)

parallel -j $CPU bwa-mem2 index {} ::: $(ls genome/*.fasta)

