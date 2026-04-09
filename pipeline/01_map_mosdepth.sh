#!/usr/bin/bash -l
#SBATCH -p short -c 24 --mem 24gb --out logs/bwamem.%a.log


CPU=${SLURM_CPUS_ON_NODE}
if [ -z $CPU ]; then
    CPU=2
fi
IFS=,
N=${SLURM_ARRAY_TASK_ID}
if [ -z $N ]; then
    N=$1
    if [ -z $N ]; then
        echo "cannot find a cmdline option or array/-a option"
        exit
    fi
fi
echo "N is $N"
SAMPLES=samples.csv
OUT=aln
DB=genome
IN=reads
RESULT=results/mosdepth

mkdir -p $OUT $RESULT
tail -n +2 $SAMPLES | sed -n ${N}p | while read ID FWD REV
do
    if [ ! -s $OUT/${ID}.bam ]; then
        module load bwa-mem2
        module load samtools
        bwa-mem2 mem -o ${SCRATCH}/${ID}.sam -t ${CPU} $DB/${ID}.masked.fasta $IN/$FWD $IN/$REV
        samtools view -OBAM -F 12 -o ${SCRATCH}/${ID}.bam ${SCRATCH}/${ID}.sam
        samtools sort -OBAM -@${CPU} -o $OUT/${ID}.bam ${SCRATCH}/${ID}.bam
    fi

    if [ ! -f $OUT/$ID.bam.bai ]; then
        samtools index $OUT/${ID}.bam
    fi
    module load mosdepth
    mosdepth -x -t $CPU $RESULT/$ID $OUT/$ID.bam
done
