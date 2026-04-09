#!/usr/bin/bash -l
#SBATCH -p short 

INDIR=download
OUTDIR=bed
for file in $(ls $INDIR/*/*.gff3.gz)
do
    name=$(basename $file .gff3.gz | perl -p -e 's/Basidiobolus_//; s/_/./')
    zgrep -P "\tgene\t" $file | cut -f1,4,5,9 | perl -p -e 's/ID=([^;]+);\S*/$1/' > $OUTDIR/${name}.bed   
done