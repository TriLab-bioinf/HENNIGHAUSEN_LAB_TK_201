#!/bin/sh

set -o errexit

source $CONFIG
mkdir -p results/8a-variant-recal-snp
threads=8

SAMPLE=$1

# load modules 
module load GATK/4.3.0.0

# STEP 1: VariantRecalibrator

# Check if sample was already processed
if [[ ! -e ./track/ALL.var_recal_snp.OK ]]; then 

    gatk --java-options "-Xmx4g -Xms4g -DGATK_STACKTRACE_ON_USER_EXCEPTION=true" VariantRecalibrator \
        --variant results/7-genotype-gvcf/all.vcf.gz \
        --reference $GENOME \
        --resource:hapmap,known=false,training=true,truth=true,prior=15.0 $HAPMAP \
        --resource:omni,known=false,training=true,truth=false,prior=12.0 $_1000G_OMNI \
        --resource:1000G,known=false,training=true,truth=false,prior=10.0 $_1000G_SNP \
        --resource:dbsnp,known=true,training=false,truth=false,prior=2.0 $DBSNP \
        -mode SNP \
        --output results/8a-variant-recal-snp/all.snp.recal.vcf \
        --tranches-file results/8a-variant-recal-snp/all.tranches \
        --rscript-file results/8a-variant-recal-snp/output.plots.R \
        --max-gaussians 8 \
        -an MQ -an QD -an MQRankSum -an ReadPosRankSum -an FS -an SOR

fi

touch ./track/ALL.var_recal_snp.OK
echo; echo Variant recalibration SNPs OK 

# STEP 2: apply_vqsr

mkdir -p results/9a-apply-recal-snp

# Check if sample was already processed
if [[ ! -e ./track/ALL.apply_vqsr_snp.OK ]]; then 

    gatk --java-options "-Xmx4g -Xms4g -DGATK_STACKTRACE_ON_USER_EXCEPTION=true" ApplyVQSR \
            --variant results/7-genotype-gvcf/all.vcf.gz \
            --recal-file results/8a-variant-recal-snp/all.snp.recal.vcf \
            --reference $GENOME \
            --tranches-file results/8a-variant-recal-snp/all.tranches \
            --mode SNP \
            --truth-sensitivity-filter-level 99.5 \
            --output results/9a-apply-recal-snp/all.snp.filtered.vcf.gz

fi

touch ./track/ALL.apply_vqsr_snp.OK
echo; echo Apply VQSR SNPs OK 

exit 0





