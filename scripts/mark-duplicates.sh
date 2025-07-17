#!/bin/sh

set -o errexit

source $CONFIG
mkdir -p results/3-mark-dup
threads=8

SAMPLE=$1

# Check if sample was already processed
if [[ -e ./track/$SAMPLE.dedup.OK ]]; then exit 0; fi 

# load modules
module load picard/2.9.2

java -Xmx60g -XX:ParallelGCThreads=${threads} -jar $PICARDJARPATH/picard.jar MarkDuplicates \
    INPUT=./results/2-mapping/$SAMPLE.bam \
    OUTPUT=./results/3-mark-dup/$SAMPLE.dedup.bam \
    METRICS_FILE=./results/3-mark-dup/${SAMPLE}.metrics.txt

java -Xmx60g -jar $PICARDJARPATH/picard.jar BuildBamIndex INPUT=./results/3-mark-dup/$SAMPLE.dedup.bam

touch ./track/$SAMPLE.dedup.OK
echo; echo Mark-duplicates $SAMPLE OK 

exit 0