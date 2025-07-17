#!/bin/sh

set -o errexit

source $CONFIG
mkdir -p results/7-genotype-gvcf
threads=8

SAMPLE=$1

# load modules 
module load GATK/4.3.0.0

# Check if sample was already processed
if [[ ! -e ./track/ALL.gen_gvcfs.OK ]]; then 

    gatk --java-options "-Xmx4g -Xms4g -DGATK_STACKTRACE_ON_USER_EXCEPTION=true" GenotypeGVCFs \
        --variant gendb://results/6-snp-db \
        --reference $GENOME \
        --intervals $WHOLE_CHROM_INTERVALS \
        -A StrandBiasBySample \
        --dbsnp $DBSNP \
        --output results/7-genotype-gvcf/all.vcf.gz

fi

touch ./track/ALL.gen_gvcfs.OK
echo; echo Loading SNP DB ALL SAMPLES OK 

exit 0





