#!/bin/bash
#
#SBATCH --mail-user=kh455@duke.edu
#SBATCH --mail-type=BEGIN,END,NONE,FAIL,REQUEUE
##SBATCH -p common                # Partition to submit to (comma separated)
#SBATCH -J splitRef1         # Job name
#SBATCH -n 1                     # Number of cores
#SBATCH -N 1                     # Ensure that all cores are on one machine
#SBATCH -t 24-00:00                # Runtime in D-HH:MM (or use minutes)
#SBATCH --mem 10000               # Memory in MB
#SBATCH -o _alignSplitRefs_%j.out # File for STDOUT (with jobid = %j) 
#SBATCH -e alignSplitRefs_%j.err       # File for STDERR (with jobid = %j)   


module load cutadapt/1.8.3-gcb01 
module load fastqc/0.11.5-fasrc01
module load java/1.7.0_60-fasrc01
module load jdk/9-gcb01

path_bbmap=/data/taylorlab/software/bbmap
export PATH=$path_bbmap:$PATH

snakemake out
#scripts/step1_splitSyncReadsMultiRef.pl 2 refs/AMA/AMA.fasta,refs/CSP/CSP.fasta out raw_reads/1 raw_reads/2 adapters/forwardPrimers.fasta adapters/reversePrimers.fasta

