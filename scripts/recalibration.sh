#!/bin/sh

set -o errexit

source $CONFIG
mkdir -p results/4-recalibration
threads=8

SAMPLE=$1

# Check if sample was already processed
if [[ -e ./track/$SAMPLE.recalib.OK ]]; then exit 0; fi 

# load modules  
module load GATK/4.3.0.0 samtools/1.17

gatk BaseRecalibrator \
    -R $GENOME \
    -I ./results/3-mark-dup/$SAMPLE.dedup.bam \
    --known-sites $DBSNP \
    -O ./results/4-recalibration/$SAMPLE.BaseRecalibrator.table

gatk ApplyBQSR \
        -R $GENOME \
        -I ./results/3-mark-dup/$SAMPLE.dedup.bam \
        --bqsr-recal-file ./results/4-recalibration/$SAMPLE.BaseRecalibrator.table \
        -O ./results/4-recalibration/$SAMPLE.recalib.bam \
        --create-output-bam-index true

# Rerun base recalibrator on recalibrated bam for evaluating recalibration efficiency
gatk BaseRecalibrator \
    -R $GENOME \
    -I ./results/4-recalibration/$SAMPLE.recalib.bam \
    --known-sites $DBSNP \
    -O ./results/4-recalibration/$SAMPLE.BaseRecalibratorAfter.table

# Generate before-after recalibratin plots
gatk AnalyzeCovariates \
    -before ./results/4-recalibration/$SAMPLE.BaseRecalibrator.table \
    -after  ./results/4-recalibration/$SAMPLE.BaseRecalibratorAfter.table \
    -plots  ./results/4-recalibration/$SAMPLE.realignedRecalPlots.pdf

touch ./track/$SAMPLE.recalib.OK
echo; echo Recalibration of $SAMPLE OK 

exit 0