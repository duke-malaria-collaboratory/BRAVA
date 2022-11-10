#!/usr/bin/perl
# Written by Joe Saelens
# Updated on 11/20/2020 to include read pair synchronization
# Updated on 10/25/2022 by Kathie Huang - converted to Python and split into 3 separate scripts

import os
import subprocess

# defining functions

def makeOutDirs(out):
    os.system('mkdir -p {}/cut'.format(out))
    os.system('mkdir -p {}/cut/1'.format(out))
    os.system('mkdir -p {}/cut/2'.format(out))
    print ("\tAdapter cut reads temporarily stored in {}/cut".format(out))

    os.system('mkdir -p {}/trim'.format(out))
    os.system('mkdir -p {}/trim/1'.format(out))
    os.system('mkdir -p {}/trim/2'.format(out))
    os.system('mkdir -p {}/trim/singleton'.format(out))
    os.system('mkdir -p {}/trim/summaries'.format(out))
    print ("\tTrimmed reads stored in /{}/trim".format(out))

def getReads(readsDir):
    if os.path.exists(readsDir):
        reads = os.listdir(readsDir)
    else:
        sys.exit("\t**Reads not found**")
    reads.sort()
    return reads

def trimReads(read1Dir, read2Dir, outDir, reads1, reads2, forward, reverse):
    size = len(reads1)
    size2 = len(reads2)
    if size != size2:
        sys.exit("\t***Paired-end read file names unequal***")

    for i in range(size):
        currRead1 = reads1[i]
        PRE = currRead1.split('_')
        prefix = PRE[0]
        print("\t Trimming {} reads".format(prefix))
        currRead2 = reads2[i]
        PRE2 = currRead2.split('_')
        prefix2 = PRE2[0]
        if prefix == prefix2:
            forwardElem = forward[0]
            reverseElem = reverse[0]
            os.system('cutadapt -g file:{0} -G file:{1} -o {2}/cut/1/{3}.1.fastq.gz -p {2}/cut/2/{3}.2.fastq.gz {4}/{5} {6}/{7}'.format(forwardElem, reverseElem, outDir, prefix, read1Dir, currRead1, read2Dir, currRead2))

            os.system('trimmomatic PE -phred33 -summary {0}/trim/summaries/{1}.summary {0}/cut/1/{1}.1.fastq.gz {0}/cut/2/{1}.2.fastq.gz {0}/trim/1/{1}.1.fastq.gz {0}/trim/singleton/{1}.1_unpaired.fq.gz {0}/trim/2/{1}.2.fastq.gz {0}/trim/singleton/{1}.2_unpaired.fq.gz LEADING:10 TRAILING:10 SLIDINGWINDOW:4:15 MINLEN:80'.format(outDir, prefix))

            os.system('rm {0}/cut/1/{1}.1.fastq.gz {0}/cut/2/{1}.2.fastq.gz'.format(outDir, prefix))

    trimReads1 = "{}/trim/1".format(outDir)
    trimReads2 = "{}/trim/2".format(outDir)
    return trimReads1, trimReads2
    
out = snakemake.params["out"]
pair1 = snakemake.params["pair1"]
pair2 = snakemake.params["pair2"]
forward = snakemake.params["forward"]
reverse = snakemake.params["rev"]

out = out.replace('\n', '')    # remove '\n' only
makeOutDirs(out);

pairedReadFiles1 = getReads(pair1)
pairedReadFiles2 = getReads(pair2)

trim1, trim2 = trimReads(pair1, pair2, out, pairedReadFiles1, pairedReadFiles2, forward, reverse)
trimReads1 = getReads(trim1)
trimReads2 = getReads(trim2)
os.system('rm -r {}/cut'.format(out))