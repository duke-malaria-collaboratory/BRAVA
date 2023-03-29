#!/bin/bash
#SBATCH --mail-type=END                 # Sends email when job is completed
#SBATCH --mail-user=<YOUR_EMAIL_HERE>   # Email of user
##SBATCH -p common                      # Partition to submit to (comma separated)
#SBATCH -J haplotype_calling_pipeline   # Job name
#SBATCH -n 1                            # Number of cores
#SBATCH -N 1                            # Ensure that all cores are on one machine
#SBATCH -t 24-00:00                     # Runtime in D-HH:MM (or use minutes)
#SBATCH --mem 10000                     # Memory in MB
#SBATCH -o logs/output_%j.out           # File for STDOUT (with jobid = %j) 
#SBATCH -e logs/output_%j.err           # File for STDERR (with jobid = %j)   

echo "$SLURM_ARRAY_TASK_ID"
snakemake --cores 12 --stats output/stats