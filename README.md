#################################
### BEFORE RUNNING
#################################

1) Start and interactive session
```
sinteractive
```

2) Edit the config.txt file 

3) Edit the input_1.txt (for read-1 files) and input_2.txt (for read-2 files) files so they point to the fastq sequencing files to be processed
   One sample per row 

4) Set the CONFIG variable pointing to the configuration file
```
CONFIG=config.txt
```

#################################
### GATK
#################################

1) FASTQC
```
swarm  --sbatch "--export=CONFIG_FILE=${CONFIG}" --logdir ./00-swarm-log --job-name step1-fastqc  -f ./scripts/WGS-pipeline-step1-fastqc.swarm
```

2) MAPPING
```
swarm  --sbatch "--export=CONFIG_FILE=${CONFIG}" --logdir ./00-swarm-log --job-name step2-mapping -f ./scripts/WGS-pipeline-step2-mapping.swarm
```

3) SPLIT CHR
```
swarm  --sbatch "--export=CONFIG_FILE=${CONFIG}" --logdir ./00-swarm-log --job-name step3-split -f ./scripts/WGS-pipeline-step3-split.swarm
```

4) DELETE DISCORDANT READS
```
swarm  --sbatch "--export=CONFIG_FILE=${CONFIG}" --logdir ./00-swarm-log --job-name step4-cleanBAM -f ./scripts/WGS-pipeline-step4-cleanBAM.swarm
```

5) MARK DUPLICATES
```
swarm  --sbatch "--export=CONFIG_FILE=${CONFIG}" --logdir ./00-swarm-log --job-name step5-MarkDuplicates -f ./scripts/WGS-pipeline-step5-MarkDuplicates.swarm
```

6) BASE RECALIBRATION
```
swarm  --sbatch "--export=CONFIG_FILE=${CONFIG}" --logdir ./00-swarm-log --job-name step6-recalibration -f ./scripts/WGS-pipeline-step6-recalibration.swarm
```

7) SPLIT RECALIBRATED BAM FILES
```
swarm  --sbatch "--export=CONFIG_FILE=${CONFIG}" --logdir ./00-swarm-log --job-name step7-split -f ./scripts/WGS-pipeline-step7-split.swarm
```

8) VARIANT CALLING
```
swarm  --sbatch "--export=CONFIG_FILE=${CONFIG}" --logdir ./00-swarm-log --job-name step8-VariantCalling -f ./scripts/WGS-pipeline-step8-VariantCalling.swarm
```

9) Filter VGF files
```
swarm  --sbatch "--export=CONFIG_FILE=${CONFIG}" --logdir ./00-swarm-log --job-name step9-filterVariants -f ./scripts/WGS-pipeline-step9-filterVariants.swarm
```


#################################
### IF YOU WANT TO CHANGE THE RUNTIME
#################################
```
newwall --jobid JOBID --time HH:MM:SS
```

#################################
### IF YOU WANT TO CANCEL YOUR JOB
#################################
```
scancel JOBID
```

# HENNIGHAUSEN_TK_201
