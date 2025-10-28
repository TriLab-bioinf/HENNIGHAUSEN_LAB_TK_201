#!/bin/bash
set -o errexit

# Auxiliary functions
check_file(){
    if [[ ! -e $1 ]]; then 
        echo; echo File $1 missing!!; echo; exit 1; 
    fi
}

# Set default configuration file to config.txt
# Otherwise you can change the name of the configuration file used to for example "my_custom_configuration_file.txt" by running this scripts like this:
# $ run_wgs_workflow.sh my_custom_configuration_file.txt

CONFIG=$1
if [[ -z $1 ]]; then
    CONFIG=${PWD}/data/config.txt
fi
check_file $CONFIG

# Load config variables
source ${CONFIG}

# Check all required files are present
for file in $METADATA $GENOME $DICT $DBSNP $HAPMAP \
            $_1000G_OMNI $_1000G_SNP $MILLS $AXIOM $BWADB_MEM2.ann $INTERVALS $WHOLE_CHROM_INTERVALS; do
    check_file $file
done

# Parse metadata file

# Sample sheet for Illumina sequencing
# sample_ID,replicate_ID,flowcell_ID,lane,library,technology,fastq_1,fastq_2
# 170220-2,A,NB500,001,1,Illumina,170220-2_R1.test.fastq.gz,170220-2_R2.test.fastq.gz

module load csvkit 2>/dev/null
samples=( $(csvcut -c sample_ID -K 14 ${METADATA})); unset samples[0]
replicate=( $(csvcut -c replicate_ID -K 14 ${METADATA})); unset replicate[0]
flowcell=( $(csvcut -c flowcell_ID -K 14 ${METADATA})); unset flowcell[0]
lane=( $(csvcut -c lane -K 14 ${METADATA})); unset lane[0]
library=( $(csvcut -c library -K 14 ${METADATA})); unset library[0]
platform=( $(csvcut -c technology -K 14 ${METADATA})); unset platform[0]
fq1=( $(csvcut -c fastq_1 -K 14 ${METADATA})); unset fq1[0]
fq2=( $(csvcut -c fastq_2 -K 14 ${METADATA})); unset fq2[0]

# Number of samples
number_samples=${#samples[@]}

# Create a (hopefully) unique prefix for the names of all jobs in this 
# particular run of the pipeline. This makes sure that runs can be
# identified unambiguously
run=$(uuidgen | tr '-' ' ' | awk '{print $1}')


# echo; echo CHECK csv parsing:; echo
# echo samples=${samples[@]}
# echo replicate=${replicate[@]}
# echo flowcell=${flowcell[@]}
# echo lane=${lane[@]}
# echo library=${library[@]}
# echo platform=${platform[@]}
# echo fq1=${fq1[@]}
# echo fq2=${fq2[@]}
# echo number_samples=$number_samples
# echo run=$run

################################################################################
#                               run the pipeline                               #
################################################################################
mkdir -p log track

echo; echo Start submitiing jobs...; echo

# STEP 1: TRIMMING. All trimming will start in parallel

declare -a trim_jobids

for idx in $(seq 1 $number_samples); do
    job_idx=$( echo $(expr ${idx} - 1))
    n=${samples[$idx]}
    echo Submitting trimming job $n using ${fq1[$idx]} ${fq2[$idx]}
    trim_jobids+=($(sbatch --cpus-per-task=8 --mem=16g --time=8:00:00 \
        --job-name=$run.trim \
        --export=CONFIG=${CONFIG} \
        --output=log/$run.trim.$n \
        scripts/trimming.sh $n ${READS}/${fq1[$idx]} ${READS}/${fq2[$idx]}))
done

# STEP 2: ALIGNING. All alignments will start in parallel 

declare -a map_jobids

for idx in $(seq 1 $number_samples); do
    job_idx=$( echo $(expr $idx - 1))
    n=${samples[$idx]}
    # @RG line
    rg_line_id="ID:${samples[$idx]}"
    rg_line_sm="SM:${samples[$idx]}"
    rg_line_pl="PL:${platform[$idx]}"
    rg_line_lb="LB:${library[$idx]}"
    
    echo Submitting mapping job for $n after trimming job ${trim_jobids[$job_idx]}
    map_jobids+=($(sbatch --cpus-per-task=8 --mem=64g --time=8:00:00 \
        --job-name=$run.align \
        --output=log/$run.align.$n \
        --export=CONFIG=${CONFIG},ID=${rg_line_id},SM=${rg_line_sm},PL=${rg_line_pl},LB=${rg_line_lb} \
        --dependency=afterok:${trim_jobids[$job_idx]} \
        scripts/align.sh $n))
done

# STEP 3: MARKING DUPLICATES. All dedup will start in parallel 

declare -a dedup_jobids

for idx in $(seq 1 $number_samples); do
    job_idx=$( echo $(expr $idx - 1))
    n=${samples[$idx]}
    echo Submitting mark duplicates job for $n after mapping job ${map_jobids[$job_idx]}
    dedup_jobids+=($(sbatch --cpus-per-task=8 --mem=64g --time=8:00:00 \
        --job-name=$run.dedup \
        --output=log/$run.mark-duplicates.$n \
        --export=CONFIG=${CONFIG} \
        --dependency=afterok:${map_jobids[$job_idx]} \
        scripts/mark-duplicates.sh $n))
done

# STEP 4: BASE RECALIBRATION.

declare -a recalib_jobids

for idx in $(seq 1 $number_samples); do
    job_idx=$( echo $(expr $idx - 1))
    n=${samples[$idx]}
    echo Submitting recalibration job for $n after mapping job ${dedup_jobids[$job_idx]}
    recalib_jobids+=($(sbatch --cpus-per-task=8 --mem=64g --time=8:00:00 \
        --job-name=$run.recal \
        --output=log/$run.recal.$n \
        --export=CONFIG=${CONFIG} \
        --dependency=afterok:${dedup_jobids[$job_idx]} \
        scripts/recalibration.sh $n))
done

# STEP 5: GCVF HAPLOTYPE CALLER

declare -a hapcaller_jobids

for idx in $(seq 1 $number_samples); do
    job_idx=$( echo $(expr $idx - 1))
    n=${samples[$idx]}
    echo Submitting variant-calling job for $n after mapping job ${recalib_jobids[$job_idx]}
    hapcaller_jobids+=($(sbatch --cpus-per-task=8 --mem=64g --time=8:00:00 \
        --job-name=$run.hapcall \
        --output=log/$run.hapcall.$n \
        --export=CONFIG=${CONFIG} \
        --dependency=afterok:${recalib_jobids[$job_idx]} \
        scripts/variant-calling.sh $n))
done

# STEP 6: UPLOAD GENOMICSDBIMPORT

declare -a snpdb_jobids

old_IFS=$IFS
IFS=":"
hapcaller_all_jobs="${hapcaller_jobids[*]}"
IFS=$old_IFS
echo Submitting DBimport load after mapping jobs ${hapcaller_jobids[@]}
snpdb_jobids+=($(sbatch --cpus-per-task=8 --mem=64g --time=8:00:00 \
    --job-name=$run.snpdb \
    --output=log/$run.snpdb \
    --export=CONFIG=${CONFIG} \
    --dependency=afterok:"${hapcaller_all_jobs}" \
    scripts/genomicsDBimport.sh))

# STEP 7: JOINED GENOTYPING USING GenomicsDB

declare -a gengvcfs_jobids

echo Submitting genotyping GVCFs after DBImport load job ${snpdb_jobids[0]}
gengvcfs_jobids+=($(sbatch --cpus-per-task=8 --mem=64g --time=8:00:00 \
    --job-name=$run.gengvcfs \
    --output=log/$run.gengvcfs \
    --export=CONFIG=${CONFIG} \
    --dependency=afterok:"${snpdb_jobids[0]}" \
    scripts/genotype_gvcfs.sh))


# STEP 8a: VARIANT RECALIBRATION AND FILTERING FOR SNPs

declare -a var_recal_snps_jobids

echo Submitting variant recallibrationn SNPs after genotype GVCFs job ${gengvcfs_jobids[0]}
var_recal_snps_jobids+=($(sbatch --cpus-per-task=8 --mem=64g --time=8:00:00 \
    --job-name=$run.varecalsnp \
    --output=log/$run.varecalsnp \
    --export=CONFIG=${CONFIG} \
    --dependency=afterok:"${gengvcfs_jobids[0]}" \
    scripts/variant_recalibrator_snp.sh))

# STEP 8b: VARIANT RECALIBRATION AND FILTERING FOR INDELs

declare -a var_recal_indel_jobids

echo Submitting variant recallibrationn INDELs after genotype GVCFs job ${gengvcfs_jobids[0]}
var_recal_indel_jobids+=($(sbatch --cpus-per-task=8 --mem=64g --time=8:00:00 \
    --job-name=$run.varecalindel \
    --output=log/$run.varecalindel \
    --export=CONFIG=${CONFIG} \
    --dependency=afterok:"${gengvcfs_jobids[0]}" \
    scripts/variant_recalibrator_indel.sh))

