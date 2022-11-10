#!/usr/bin/perl
# Written by Joe Saelens
# Updated on 11/20/2020 to include read pair synchronization
# Updated on 10/25/2022 by Kathie Huang - converted to Python and split into 3 separate scripts

import os
import subprocess

# defining functions

def makeOutDirs(out):
    i = 0
    os.system('mkdir -p {}'.format(out))
    os.system('mkdir -p {}/fastqc_in'.format(out))
    print ("\tFastqc .html files of input reads stored in {}/fastqc_in".format(out))

def getReads(readsDir):
    if os.path.exists(readsDir):
        reads = os.listdir(readsDir)
    else:
        sys.exit("\t**Reads not found**")
    reads.sort()
    return reads

def QCreads(read1Dir, read2Dir, outDir, reads1, reads2, formatting):
    size = len(reads1)
    size2 = len(reads2)
    if size != size2:
        sys.exit("Paired end read files unequal")
    for i in range(size):
        os.system('fastqc -o {0} -f {1} {2}/{3} {4}/{5}'.format(outDir, formatting, read1Dir, reads1[i], read2Dir, reads2[i]))

# retrieving configurable paths
out = snakemake.params["out"]
pair1 = snakemake.params["pair1"]
pair2 = snakemake.params["pair2"]

out = out.replace('\n', '') # remove '\n' only

makeOutDirs(out);

pairedReadFiles1 = getReads(pair1)
pairedReadFiles2 = getReads(pair2)

# call fastqc
QCreads(pair1, pair2, out + "/fastqc_in", pairedReadFiles1, pairedReadFiles2, "fastq")