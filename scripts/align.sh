#!/bin/bash

set -o errexit

source $CONFIG
mkdir -p results/2-mapping
threads=8

SAMPLE=$1

# Check if sample was already processed
if [[ -e ./track/$SAMPLE.align.OK ]]; then exit 0; fi 

module load bwa-mem2/2.2.1 samtools/1.21 2>/dev/null
        
bwa-mem2 mem \
    -t $threads \
    $BWADB_MEM2 ./results/1-trimming/$SAMPLE.P.R1.fq.gz ./results/1-trimming/$SAMPLE.P.R2.fq.gz | \
samtools view -hb - | \
samtools addreplacerg -r $ID -r $SM -r $PL -r $LB  - | \
samtools sort -@ $threads \
    -O BAM -T ./results/2-mapping/$SAMPLE.tmp \
    --write-index \
    -o ./results/2-mapping/$SAMPLE.bam -

samtools index -@ $threads ./results/2-mapping/$SAMPLE.bam

touch ./track/$SAMPLE.align.OK
echo; echo Aligning $SAMPLE OK 

exit 0
