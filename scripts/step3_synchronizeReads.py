#!/usr/bin/perl
# Written by Joe Saelens
# Updated on 11/20/2020 to include read pair synchronization
# Updated on 10/25/2022 by Kathie Huang - converted to Python and split into 3 separate scripts

import os
import subprocess
import re

# defining functions

def getReference(numRef, refPaths):
    i = 0
    print(numRef)
    print(refPaths)
    
    refs = refPaths.split(',')
    refNames = [None] * len(refs)
    referenceFastas = [None] * len(refs)
    for ref in refs:
        referenceFastas[i] = ref
        # remove leading and trailing whitespace
        referenceFastas[i] = referenceFastas[i].strip()
        if not os.path.exists(referenceFastas[i]):
            sys.exit("**Reference file does not exist**")
        refSplit = referenceFastas[i].split('/')
        x = len(refSplit)
        nameIndex = x - 1
        refFasta = refSplit[nameIndex].split('\.')
        refNames[i] = refFasta[0]
        i += 1
    return referenceFastas, refNames

def makeOutDirs(out):
    os.system('mkdir -p {}/fastq/'.format(out))
    print("\tSynchronized paired end fastq files stored in /{}/fastq".format(out))
    os.system('mkdir -p {}/results'.format(out))
    print("\tLog files stored in /{}/results".format(out))

def getReads(readsDir):
    if os.path.exists(readsDir):
        reads = os.listdir(readsDir)
    else:
        sys.exit("\t**Reads not found**")
    reads.sort()
    return reads

def splitReads(refSeqs, refNames, out, read1, read2, read1Dir, read2Dir, target):
    size = len(read1)
    refs = ",".join(refSeqs)
    print("refs")
    print(refs)
    sampleName = [None] * size
    for i in range(size):
        fileSplit = re.split('\.', read1[i])
        sampleName[i] = fileSplit[0]
        print("\tAligning {0} to {1}".format(sampleName[i], refs))
        file2Split = re.split('\.', read2[i])
        sampleName2 = file2Split[0]

        if sampleName[i] == sampleName2:
            os.system('bbsplit.sh -Xmx15g in={0}/{1} in2={2}/{3} ref={4} basename={5}_%_#.fastq >& {6}/results/{5}.txt overwrite=true nodisk path={6}'.format(read1Dir, read1[i], read2Dir, read2[i], refs, sampleName[i], out))

    for elem in refNames:
        os.system("mv *_{0}* {1}/fastq".format(target, out))
        os.system("gzip {0}/fastq/*.fastq".format(out))
        if not os.path.exists(out + "/fastq/1"):
            os.system("mkdir {}/fastq/1".format(out))
        if not os.path.exists(out + "/fastq/2"):
            os.system("mkdir {}/fastq/2".format(out))
        os.system("mv {0}/fastq/*_1.fastq.gz {0}/fastq/1".format(out))
        os.system("mv {0}/fastq/*_2.fastq.gz {0}/fastq/2".format(out))
    print("Synchronizing paired end reads...\n")
    for elem in refNames:
        unsyncReads1 = getReads("{}/fastq/1".format(out))
        unsyncReads2 = getReads("{}/fastq/2".format(out))
        syncReads("{}/fastq/1".format(out), "{}/fastq/2".format(out), unsyncReads1, unsyncReads2)

def syncReads(read1Dir, read2Dir, reads1, reads2):
    size = len(reads1)
    size2 = len(reads2)
    if size != size2:
        sys.exit("Paired end read files unequal")
    sampleName = [None] * size

    for i in range(size):
        fileSplit = re.split('\.', reads1[i])
        sampleName[i] = fileSplit[0]
        print("\tSynchronizing paired read ends for {}".format(sampleName[i]))
        file2Split = re.split('\.', reads2[i])
        sampleName2 = file2Split[0]

        os.system("repair.sh in={0}/{1} in2={2}/{3} out={0}/{4}.paired.fastq.gz out2={2}/{5}.paired.fastq.gz outs=singletons.fq repair".format(read1Dir, reads1[i], read2Dir, reads2[i], sampleName[i], sampleName2))
        os.system("rm singletons.fq {0}/{1} {2}/{3}".format(read1Dir, reads1[i], read2Dir, reads2[i]))
        os.system("mv {0}/{1}.paired.fastq.gz {0}/{2}".format(read1Dir, sampleName[i], reads1[i]))
        os.system("mv {0}/{1}.paired.fastq.gz {0}/{2}".format(read2Dir, sampleName2, reads2[i]))

numRef = 1
refs = snakemake.params["refs"]
out = snakemake.params["out"]
target = snakemake.params["target"]

refSeqs, refNames = getReference(numRef, refs)

out = out.replace('\n', '')    # remove '\n' only
makeOutDirs(out);

trim1 = snakemake.params["trimmed"] + "/1"
trim2 = snakemake.params["trimmed"] + "/2"

trimReads1 = getReads(trim1)
trimReads2 = getReads(trim2)

splitReads(refSeqs, refNames, out, trimReads1, trimReads2, trim1, trim2, target)

os.system('mkdir {}'.format(snakemake.params["all_samples"]))
os.system('mv {0}/*.fastq.gz {1} && mv {2}/*.fastq.gz {1}'.format(snakemake.params["forward_samples"], snakemake.params["all_samples"], snakemake.params["reverse_samples"]))
os.system('rm -r {0} && rm -r {1}'.format(snakemake.params["forward_samples"], snakemake.params["reverse_samples"]))