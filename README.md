# BASH-based WES/WGS SNP Identification Pipeline Quick Start (Biowulf Cluster)

How to run the BASH-based Whole Exome/Genome Sequencing SNP identification pipeline on NIH's Biowulf cluster.

## Prerequisites

- Access to Biowulf cluster
- FASTQ files from Illumina sequencing
- Basic knowledge of SLURM job submission

## Quick Setup

### 1. Download and Setup

```bash
# Navigate to your data directory
cd /data/$USER

# Clone/download the pipeline
git clone https://github.com/TriLab-bioinf/HENNIGHAUSEN_LAB_TK_201.git wgs-pipeline
cd wgs-pipeline

# Make scripts executable
chmod +x run_wgs_workflow.sh scripts/*.sh
```

### 2. Configure Sample Sheet

Edit `data/samplesheet.csv` with your samples:

```csv
sample_ID,replicate_ID,flowcell_ID,lane,library,technology,fastq_1,fastq_2
SAMPLE1,A,FC001,001,1,Illumina,SAMPLE1_R1.fastq.gz,SAMPLE1_R2.fastq.gz
SAMPLE2,A,FC001,002,1,Illumina,SAMPLE2_R1.fastq.gz,SAMPLE2_R2.fastq.gz
```

### 3. Update Configuration

Edit `data/config.txt` - **only change these lines**:

```bash
# Update this path to your FASTQ files location
READS=/data/$USER/path/to/your/fastq/files

# For first run (create new database)
DB_ACTION=create

# For adding samples to existing results
# DB_ACTION=update
```

**Note**: All reference files and modules are pre-configured for Biowulf.

## Run Pipeline

### First Time Run

```bash
# Submit the pipeline
./run_wgs_workflow.sh

# Monitor jobs
squeue -u $USER
```

### Adding New Samples

1. Add new samples to `data/samplesheet.csv`
2. Change `DB_ACTION=update` in `data/config.txt`
3. Run: `./run_wgs_workflow.sh`

## Monitor Progress

```bash
# Check running jobs
squeue -u $USER

# View log files
ls log/
tail -f log/<run_id>.trim.SAMPLE1

# Check completed steps
ls track/
```

## Key Directories and Output Files

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

## Troubleshooting

### Common Issues

1. **Job failures**: Check logs in `log/` directory
2. **File not found**: Verify `READS` path in `config.txt`
3. **Permission errors**: Ensure scripts are executable (`chmod +x`)

### Restart Failed Jobs

The pipeline automatically resumes from the last completed step. To force rerun:

```bash
# Rerun specific sample step
rm track/SAMPLE1.trim.OK

# Rerun all steps for a sample
rm track/SAMPLE1.*.OK

# Start completely fresh
rm track/*.OK
```

---

For detailed documentation, see [DOCUMENTATION.md](./DOCUMENTATION.md)

## Quick Commands Reference

```bash
# Check job status
squeue -u $USER

# Cancel all jobs
scancel -u $USER

# Check disk usage
du -skh .

# Check sample sheet
module load csvkit
csvstat -K 14 data/samplesheet.csv

```
