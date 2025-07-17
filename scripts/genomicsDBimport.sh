#!/bin/sh

set -o errexit

source $CONFIG
threads=8

declare -a joined_inputs

# load modules 
module load GATK/4.3.0.0

# Check if samples were already processed
if [[ ! -e ./track/ALL.snpdb.OK || $DB_ACTION == "update" ]]; then 

    # Collect individual gvcf files for input
    for sample in $(/bin/ls ./results/5-variant-calls/*.g.vcf.gz); do
    joined_inputs+="--variant "
    joined_inputs+="$sample "
    done

    # Create new SNP-DB or update existent SNP-DB
    if [[ $DB_ACTION == "create" ]]; then
        DB="--genomicsdb-workspace-path results/6-snp-db"
    elif [[ $DB_ACTION == "update" ]]; then
        DB="--genomicsdb-update-workspace-path results/6-snp-db"
    else    
        echo; echo ERROR, DB_ACTION should be either 'create' or 'update'. Current value = "${DB_ACTION}" ; echo
        exit 1
    fi


    gatk --java-options "-Xmx4g -Xms4g -DGATK_STACKTRACE_ON_USER_EXCEPTION=true" GenomicsDBImport ${joined_inputs} --reader-threads 1 \
        --genomicsdb-shared-posixfs-optimizations true --merge-contigs-into-num-partitions 25 \
        -L $WHOLE_CHROM_INTERVALS $DB
        
    touch ./track/ALL.snpdb.OK
    echo; echo Loading SNP DB ALL SAMPLES OK 

    # Make SNPDB backup

    mkdir -p results/6-snp-db-bkp            
    tar -czf results/6-snp-db-bkp/genomeDB_$(echo $$).tar.gz 6-snp-db
    touch results/6-snp-db/genomeDB_id.$$

fi

exit 0
