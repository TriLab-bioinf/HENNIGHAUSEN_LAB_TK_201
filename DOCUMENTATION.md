# Whole Exome/Genome Sequencing (WES/WGS) SNP Identification Pipeline (Human genome)

This repository contains a comprehensive bash-based pipeline for processing whole genome sequencing (WGS) data using GATK best practices. The pipeline performs quality control, alignment, variant calling, and variant quality score recalibration (VQSR).

## Table of Contents

- [Overview](#overview)
- [Prerequisites](#prerequisites)
- [Pipeline Architecture](#pipeline-architecture)
- [Traking completed jobs](#tracking-of-job-completions)
- [Installation and Setup](#installation-and-setup)
- [Configuration](#configuration)
- [Running the Pipeline](#running-the-pipeline)
- [Output Files](#output-files)
- [Troubleshooting](#troubleshooting)
- [Pipeline Steps](#pipeline-steps)

## Overview

This pipeline implements the GATK Best Practices workflow for germline variant discovery in whole genome sequencing data. It includes the following major steps:

1. **Quality Control & Trimming** - Adapter removal and quality trimming
2. **Alignment** - BWA-MEM2 alignment to reference genome
3. **Mark Duplicates** - Remove PCR duplicates using GATK MarkDuplicates
4. **Base Quality Score Recalibration (BQSR)** - Correct systematic errors in base quality scores
5. **Variant Calling** - Generate GVCF files using GATK HaplotypeCaller
6. **Joint Genotyping** - Consolidate variants across all samples
7. **Variant Quality Score Recalibration (VQSR)** - Filter variants based on quality metrics

## Prerequisites

### Required Software

- **SLURM** - Workload manager for job scheduling
- **BWA-MEM2** - Fast and accurate aligner
- **GATK** (version 4.x) - Genome Analysis Toolkit
- **SAMtools** - Utilities for manipulating alignments
- **csvkit** - Command-line tools for working with CSV files
- **trimmomatic** or similar trimming software

### Required Reference Files

The pipeline requires several reference files and databases:

- Human reference genome (GRCh38/hg38)
- BWA-MEM2 index files
- GATK bundle files including:
  - dbSNP database
  - HapMap variants
  - 1000 Genomes variants
  - Mills and 1000G gold standard indels
  - Axiom Exome variants

### System Requirements

- Linux/Unix environment with SLURM scheduler
- Sufficient storage space for intermediate and output files
- Adequate compute resources (recommended: 8+ CPUs, 64GB+ RAM per job)

## Pipeline Architecture

The pipeline is designed to run on SLURM clusters and submits sbatch jobs with dependencies to ensure proper execution order:

```
Trimming (parallel) → Alignment (parallel) → Mark Duplicates (parallel) → 
Base Recalibration (parallel) → Variant Calling (parallel) → 
GenomicsDB Import → Joint Genotyping → VQSR (SNPs & INDELs in parallel)
```

## Job Completion Tracking

The pipeline implements an intelligent tracking system to monitor job completion status and enable pipeline restarts without reprocessing completed steps.

### How It Works

The pipeline creates tracking files in the `track/` directory using the following naming convention:

```
track/<sample_ID>.<step_name>.OK    # For per-sample steps
track/ALL.<step_name>.OK            # For multi-sample steps
```

### Tracking File Examples

```bash
track/SAMPLE1.trim.OK               # Sample1 trimming completed
track/SAMPLE1.align.OK              # Sample1 alignment completed
track/SAMPLE2.dedup.OK              # Sample2 duplicate marking completed
track/ALL.snpdb.OK                  # GenomicsDB import completed
track/ALL.gengvcfs.OK               # Joint genotyping completed
```

### Benefits of the Tracking System

1. **Resume Capability**: If the pipeline fails or is interrupted, it can resume from where it left off
2. **Incremental Processing**: When adding new samples, only new samples are processed while existing results are preserved
3. **Resource Efficiency**: Avoids recomputing expensive steps that have already completed successfully
4. **Manual Control**: Users can selectively rerun specific steps by removing tracking files

### Manual Tracking Management

To **rerun a specific step**, delete its corresponding tracking file:

```bash
# Rerun trimming for SAMPLE1
rm track/SAMPLE1.trim.OK

# Rerun joint genotyping for all samples
rm track/ALL.gengvcfs.OK

# Rerun all steps for SAMPLE2
rm track/SAMPLE2.*.OK
```

To **start completely fresh**, remove all tracking files:

```bash
# Remove all tracking files
rm track/*.OK
```

### Important Notes

- Tracking files are created only upon **successful completion** of each step
- The pipeline checks for tracking files before submitting jobs
- Dependencies between steps are still respected even with tracking files present
- Tracking files are small text files and safe to manually manage 

## Installation and Setup

1. **Clone or download the pipeline**:
   ```bash
   git clone <repository-url>
   cd HENNIGHAUSEN_LAB_TK_201
   ```

2. **Make scripts executable**:
   ```bash
   chmod +x run_wgs_workflow.sh
   chmod +x scripts/*.sh
   ```

3. **Required modules** (Available already in Biowulf):
   ```bash
   fastp/0.24.0
   GATK/4.3.0.0
   picard/2.9.2
   bwa-mem2/2.2.1
   samtools/1.21
   csvkit/2.1.0
   ```

## Configuration

### 1. Sample Sheet Configuration

Edit `data/samplesheet.csv` to include your samples:

```csv
sample_ID,replicate_ID,flowcell_ID,lane,library,technology,fastq_1,fastq_2
SAMPLE1,A,FC001,001,1,Illumina,SAMPLE1_R1.fastq.gz,SAMPLE1_R2.fastq.gz
SAMPLE2,A,FC001,002,1,Illumina,SAMPLE2_R1.fastq.gz,SAMPLE2_R2.fastq.gz
```

**Column descriptions**:
- `sample_ID`: Unique identifier for each sample
- `replicate_ID`: Replicate identifier (typically A, B, C, etc.)
- `flowcell_ID`: Sequencing flowcell identifier
- `lane`: Sequencing lane number
- `library`: Library preparation identifier
- `technology`: Sequencing platform (typically "Illumina")
- `fastq_1`: Forward read FASTQ filename
- `fastq_2`: Reverse read FASTQ filename

### 2. Configuration File Setup

Edit `data/config.txt` to specify file paths and parameters:

#### Essential Parameters to Update:

```bash
# Sample metadata
METADATA=./data/samplesheet.csv

# Reference files (update paths for your system)
GENOME=/path/to/Homo_sapiens_assembly38.fasta
DICT=/path/to/Homo_sapiens_assembly38.dict
DBSNP=/path/to/dbsnp_146.hg38.vcf.gz
HAPMAP=/path/to/hapmap_3.3.hg38.vcf.gz
_1000G_OMNI=/path/to/1000G_omni2.5.hg38.vcf.gz
_1000G_SNP=/path/to/1000G_phase1.snps.high_confidence.hg38.vcf.gz
MILLS=/path/to/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz
AXIOM=/path/to/Axiom_Exome_Plus.genotypes.all_populations.poly.hg38.vcf.gz

# BWA-MEM2 index
BWADB_MEM2=/path/to/bwa_mem2_index/Homo_sapiens_assembly38.fasta

# Interval files
INTERVALS=data/intervals/exons_intervals.merged.bed
WHOLE_CHROM_INTERVALS=data/intervals/whole_chr_intervals.bed

# Adapters for trimming
ADAPTERS=data/adapters/adapters.fa

# Input reads directory
READS=/path/to/your/fastq/files
```

### 3. Interval Files

Ensure interval files are properly configured:

- **`data/intervals/exons_intervals.merged.bed`**: Target intervals for variant calling (typically exonic regions for WES, whole genome intervals for WGS)
- **`data/intervals/whole_chr_intervals.bed`**: Chromosome-level intervals for GenomicsDB import

### 4. Adapter File

Ensure `data/adapters/adapters.fa` contains the appropriate sequencing adapters for your platform.

### 5. Database Action Configuration

The `DB_ACTION` parameter controls how the GenomicsDB is handled during the pipeline execution. This is crucial for managing incremental sample processing:

#### For First-Time Pipeline Execution

When running the pipeline for the first time with a new set of samples:

```bash
# In data/config.txt
DB_ACTION=create
```

**What this does**:
- Creates a new GenomicsDB from scratch
- The output directory must be empty or non-existent
- All samples in the current samplesheet will be included in the database
- **Important**: If the GenomicsDB directory already exists, the job will fail

#### For Adding New Samples to Existing Data

When you want to add new samples to an existing GenomicsDB (incremental processing):

```bash
# In data/config.txt
DB_ACTION=update
```

**What this does**:
- Updates an existing GenomicsDB with new samples
- Preserves previously processed samples
- Only processes the new samples listed in the current samplesheet
- The GenomicsDB directory must already exist from a previous run

#### Important Considerations

1. **First Run**: Always use `DB_ACTION=create` for your initial pipeline execution
2. **Adding Samples**: Use `DB_ACTION=update` when adding new samples to existing results
3. **Directory Management**: 
   - For `create`: Ensure the GenomicsDB output directory is empty/non-existent
   - For `update`: Ensure the GenomicsDB directory exists from previous runs
4. **Sample Management**: The pipeline will process all samples listed in the current `samplesheet.csv`

#### Example Workflow

**Initial run with 10 samples**:
```bash
# Set in config.txt
DB_ACTION=create

# Run pipeline
./run_wgs_workflow.sh
```

**Later, adding 5 more samples**:
```bash
# Add new samples to samplesheet.csv (can include old samples or just new ones)
# Set in config.txt
DB_ACTION=update

# Run pipeline again
./run_wgs_workflow.sh
```

### 6. Additional Configuration Options

#### Duplicate Read Handling

Control how duplicate reads are handled during processing:

```bash
# In data/config.txt
REMOVE_DUPLICATED_READS=False  # Mark duplicates but keep them
REMOVE_DUPLICATED_READS=True   # Remove duplicate reads entirely
```

**Recommendation**: Set to `False` to preserve duplicate reads for downstream analysis while still marking them.

## Running the Pipeline

### Basic Usage

Run the pipeline with default configuration:

```bash
./run_wgs_workflow.sh
```

### Using Custom Configuration

Run with a custom configuration file:

```bash
./run_wgs_workflow.sh /path/to/custom_config.txt
```

### Pre-run Checklist

Before running the pipeline, verify:

1. ✅ All reference files exist and are accessible
2. ✅ FASTQ files are in the specified `READS` directory
3. ✅ Sample sheet is correctly formatted
4. ✅ All required modules are loaded
5. ✅ Sufficient disk space is available
6. ✅ SLURM is operational and accessible

### Monitoring Jobs

Monitor job progress using SLURM commands:

```bash
# Check job queue
squeue -u $USER

# Check job details
scontrol show job <job_id>

# Monitor log files
tail -f log/<run_id>.*
```

## Output Files

The pipeline generates several types of output files:

### Directory Structure

```
├── log/                    # SLURM job logs
│   ├── <run_id>.trim.*     # Trimming logs
│   ├── <run_id>.align.*    # Alignment logs
│   ├── <run_id>.dedup.*    # Duplicate marking logs
│   ├── <run_id>.recal.*    # Recalibration logs
│   └── <run_id>.hapcall.*  # Variant calling logs
├── track/                  # Intermediate tracking files
└── results/                # Final output files (created by scripts)
    ├── trimmed/           # Quality-trimmed FASTQ files
    ├── aligned/           # Aligned BAM files
    ├── dedup/             # Deduplicated BAM files
    ├── recal/             # Recalibrated BAM files
    ├── gvcf/              # Individual GVCF files
    └── vcf/               # Final VCF files
```

### Key Output Files

- **Final VCF**: Joint-called, VQSR-filtered variant calls
- **Individual GVCFs**: Per-sample genomic variant calls
- **Recalibrated BAMs**: Analysis-ready alignment files
- **QC Reports**: Quality control metrics and statistics

## Troubleshooting

### Common Issues

1. **File Not Found Errors**:
   - Verify all paths in `config.txt` are correct and accessible
   - Check that FASTQ files exist in the `READS` directory
   - Ensure reference files are properly indexed

2. **Job Failures**:
   - Check SLURM logs in the `log/` directory
   - Verify sufficient resources are allocated
   - Check module loading in job scripts

3. **Memory Issues**:
   - Increase memory allocation in job submissions
   - Consider splitting large intervals for parallel processing

4. **Disk Space**:
   - Monitor disk usage throughout the pipeline
   - Clean up intermediate files if necessary
   - Ensure output directories have sufficient space

### Debugging Steps

1. **Check configuration**:
   ```bash
   # Verify config file syntax
   bash -n data/config.txt
   
   # Test file accessibility
   ls -la $(grep "^GENOME=" data/config.txt | cut -d'=' -f2)
   ```

2. **Validate sample sheet**:
   ```bash
   # Check CSV format
   csvstat -K 14 data/samplesheet.csv
   
   # Verify FASTQ files exist
   while IFS=, read -r sample rep fc lane lib tech fq1 fq2; do
       [[ "$sample" == "sample_ID" ]] && continue
       ls -la "$READS/$fq1" "$READS/$fq2"
   done < data/samplesheet.csv
   ```

3. **Test individual scripts**:
   ```bash
   # Test trimming script on one sample
   sbatch scripts/trimming.sh SAMPLE1 /path/to/R1.fastq.gz /path/to/R2.fastq.gz
   ```

## Pipeline Steps

### Step 1: Quality Control and Trimming
- Removes sequencing adapters
- Trims low-quality bases
- Filters short reads

### Step 2: Alignment
- Aligns reads to reference genome using BWA-MEM2
- Adds read group information
- Converts to sorted BAM format

### Step 3: Mark Duplicates
- Identifies and marks PCR duplicates
- Generates duplicate metrics

### Step 4: Base Quality Score Recalibration
- Builds recalibration model using known variant sites
- Applies recalibration to base quality scores

### Step 5: Variant Calling
- Calls variants using GATK HaplotypeCaller in GVCF mode
- Generates per-sample genomic VCF files

### Step 6: GenomicsDB Import
- Consolidates GVCF files into GenomicsDB format
- Prepares for joint genotyping

### Step 7: Joint Genotyping
- Performs joint genotyping across all samples
- Generates raw VCF with all variants

### Step 8: Variant Quality Score Recalibration
- Builds recalibration models for SNPs and INDELs separately
- Applies quality filters to variants
- Generates final filtered VCF

---

## Support

For questions or issues:

1. Check the troubleshooting section above
2. Review SLURM job logs in the `log/` directory
3. Consult GATK documentation for tool-specific issues
4. Contact your system administrator for cluster-related problems

## License

This pipeline is provided as-is for research purposes. Please cite appropriate tools and references when using this pipeline in publications.
