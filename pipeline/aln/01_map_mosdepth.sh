#!/bin/bash -l
#SBATCH -p short -c 24 --mem 24gb --out logs/bwamem.%a.log

: "${SCRATCH:?SCRATCH environment variable is not set}"

CPU=${SLURM_CPUS_ON_NODE}
if [ -z "$CPU" ]; then
    CPU=2
fi
IFS=,
N=${SLURM_ARRAY_TASK_ID}
if [ -z "$N" ]; then
    N=$1
    if [ -z "$N" ]; then
        echo "cannot find a cmdline option or array/-a option"
        exit 1
    fi
fi
echo "N is $N"
SAMPLES=samples.csv
OUT=aln
DB=genome
IN=reads
BEDIN=bed
RESULT=results/mosdepth

mkdir -p $OUT $RESULT
tail -n +2 $SAMPLES | sed -n ${N}p | while read ID FWD REV
do
    if [ ! -s $OUT/${ID}.bam ]; then
        if [[ $FWD == *pacbio* ]]; then
            module load minimap2
            module load samtools
            minimap2 -t $CPU -ax map-pb $DB/$ID.masked.fasta $IN/$FWD | samtools view -OBAM -F 12 -o $SCRATCH/${ID}.bam -
            samtools sort -OBAM -@${CPU} -o $OUT/${ID}.bam $SCRATCH/${ID}.bam
            echo "pacbiodata"
        else
            module load bwa-mem2
            module load samtools
            bwa-mem2 mem -o ${SCRATCH}/${ID}.sam -t ${CPU} $DB/${ID}.masked.fasta $IN/$FWD $IN/$REV
            samtools view -OBAM -F 12 -o ${SCRATCH}/${ID}.bam ${SCRATCH}/${ID}.sam
            samtools sort -OBAM -@${CPU} -o $OUT/${ID}.bam ${SCRATCH}/${ID}.bam
        fi
    fi

    if [ ! -f $OUT/$ID.bam.bai ]; then
        samtools index $OUT/${ID}.bam
    fi
    module load mosdepth
    BED=$BEDIN/$ID.bed
    export MOSDEPTH_Q0=NO_COVERAGE
    export MOSDEPTH_Q1=LOW_COVERAGE
    export MOSDEPTH_Q2=CALLABLE
    export MOSDEPTH_Q3=HIGH_COVERAGE
    export MOSDEPTH_Q4=VERY_HIGH_COVERAGE
    mosdepth --quantize 0:1:4:100:200: -x -t $CPU -b $BED $RESULT/$ID $OUT/$ID.bam
done
