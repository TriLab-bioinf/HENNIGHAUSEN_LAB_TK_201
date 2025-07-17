#!/bin/sh

set -o errexit

source $CONFIG
mkdir -p results/5-variant-calls
threads=8

SAMPLE=$1

# load modules 
module load GATK/4.3.0.0

# Check if sample was already processed
if [[ ! -e ./track/$SAMPLE.hapcall.OK ]]; then 

    gatk HaplotypeCaller \
        -R $GENOME \
        -I ./results/4-recalibration/$SAMPLE.recalib.bam \
        -O ./results/5-variant-calls/$SAMPLE.hapcall.g.vcf.gz \
        -ERC GVCF \
        --standard-min-confidence-threshold-for-calling 30 \
        --native-pair-hmm-threads 4 \
        -dont-use-soft-clipped-bases

    touch ./track/$SAMPLE.hapcall.OK
    echo; echo Haplotype caller on $SAMPLE OK 

fi

exit 0
