#!/bin/bash

set -o errexit

source $CONFIG
mkdir -p results/1-trimming

module load fastp/0.24.0 2>/dev/null

SAMPLE=$1        
FQ1=$2
FQ2=$3

# Check if sample was already processed
if [[ -e ./track/$SAMPLE.trim.OK ]]; then exit 0; fi 

fastp -i $FQ1 --in2 $FQ2 \
    --out1 ./results/1-trimming/$SAMPLE.P.R1.fq.gz --out2 ./results/1-trimming/$SAMPLE.P.R2.fq.gz \
    --unpaired1 ./results/1-trimming/$SAMPLE.UP.R1.fq.gz \
    --unpaired2 ./results/1-trimming/$SAMPLE.UP.R2.fq.gz \
    --failed_out ./results/1-trimming/$SAMPLE.failed.fq.gz \
    --html ./results/1-trimming/$SAMPLE.html \
    --json ./results/1-trimming/$SAMPLE.json \
    --trim_poly_x --poly_x_min_len 10 \
    --qualified_quality_phred 20 \
    --length_required 50 \
    --adapter_fasta $ADAPTERS \
    --detect_adapter_for_pe \
    --cut_right --cut_right_mean_quality 20 \
    --cut_right_window_size 5

touch ./track/$SAMPLE.trim.OK
echo; echo Mapping $SAMPLE OK 

exit 0
