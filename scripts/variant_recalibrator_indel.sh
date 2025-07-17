#!/bin/sh

set -o errexit

source $CONFIG
mkdir -p results/8b-variant-recal-indel
threads=8

SAMPLE=$1

# load modules 
module load GATK/4.3.0.0

# STEP 1: VariantRecalibrator

# Check if sample was already processed
if [[ ! -e ./track/ALL.var_recal_indel.OK ]]; then 

    gatk --java-options "-Xmx4g -Xms4g -DGATK_STACKTRACE_ON_USER_EXCEPTION=true" VariantRecalibrator \
        --variant results/7-genotype-gvcf/all.vcf.gz \
        --reference $GENOME \
        --resource:mills,known=false,training=true,truth=true,prior=12.0 $MILLS \
        --resource:axiom,known=false,training=true,truth=false,prior=10.0 $AXIOM \
        --resource:dbsnp,known=true,training=false,truth=false,prior=2.0 $DBSNP \
        -mode INDEL \
        --output results/8b-variant-recal-indel/all.indel.recal.vcf \
        --tranches-file results/8b-variant-recal-indel/all.indel.tranches \
        --rscript-file results/8b-variant-recal-indel/output.indel.plots.R \
        --max-gaussians 8 \
        -an DP -an QD -an MQRankSum -an ReadPosRankSum -an FS -an SOR

fi

touch ./track/ALL.var_recal_indel.OK
echo; echo Variant recalibration INDELs OK 

# STEP 2: apply_vqsr

mkdir -p results/9b-apply-recal-indel

# Check if sample was already processed
if [[ ! -e ./track/ALL.apply_vqsr_indel.OK ]]; then 

    gatk --java-options "-Xmx4g -Xms4g -DGATK_STACKTRACE_ON_USER_EXCEPTION=true" ApplyVQSR \
            --variant results/7-genotype-gvcf/all.vcf.gz \
            --recal-file results/8b-variant-recal-indel/all.indel.recal.vcf \
            --reference $GENOME \
            --tranches-file results/8b-variant-recal-indel/all.indel.tranches \
            -mode INDEL \
            --truth-sensitivity-filter-level 99.5 \
            --output results/9b-apply-recal-indel/all.indel.filtered.vcf.gz

fi

touch ./track/ALL.apply_vqsr_indel.OK
echo; echo Apply VQSR INDELs OK 

exit 0





