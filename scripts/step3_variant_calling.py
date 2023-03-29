#!/usr/bin/env python

import os
import subprocess
import argparse
import glob
import sys
import re

"""
This script is called as a single sample submission or Slurm array job, and takes in
the directory of the reference sequence,  forward read, reverse read, minimum base 
quality, and the output directory. If users do not submit a minumum base quality 
used for filtering variants, a phred score of 20 will be applied.

Sorted BAM and VCF files are output into the specified directory.

*************************************************************************************
Forward and reverse reads must be stored in separate directories.
*************************************************************************************

command: python variantArray.py -ref <ref_dir> -p1 <pair1> -p2 <pair2> -bq <baseQual> -o <out_dir>
"""

#Functions--------------------------------------------------------------------------
def checkVersion():
	if sys.version_info[0] < 3:
		raise Exception("Must be using Python 3")
		print("\nType: module load python\n")

def makeOutDirs(topOut):
	sam = os.path.join(topOut, "sam")
	bam = os.path.join(topOut, "bam")
	sortBam = os.path.join(topOut, "sort_bam")
	vcf = os.path.join(topOut, "vcf")
	if not os.path.exists(sam):
		os.makedirs(sam)
	if not os.path.exists(bam):
		os.makedirs(bam)
	if not os.path.exists(sortBam):
		os.makedirs(sortBam)
	if not os.path.exists(vcf):
		os.makedirs(vcf)

def getReads(readsDir):
    if os.path.exists(readsDir):
        reads = os.listdir(readsDir)
    else:
        sys.exit("\t**Reads not found**")
    reads.sort()
    return reads

def getRef(refDir):
	i = 0
	for f in os.listdir(refDir):
		if f.endswith((".fasta", ".fa")):
			fasta = os.path.join(refDir, f)
			i+=1
		elif f.endswith((".amb", ".ann", ".bwt", ".fai", ".pac", ".sa")):
			i+=1
		
	if i < 7:
		print("\n***Indexing reference sequence***\n")
		os.system("bwa index {}".format(fasta))
		os.system("samtools faidx {}".format(fasta))
	else:
		print("\n***Reference sequence already indexed***\n")
	return(fasta)

def align(refSeq, read1, read2, read1Dir, read2Dir, out):
	size = len(read1)
	sampleName = [None] * size
	for i in range(size):
		fileSplit = re.split('\.', read1[i])
		sampleName[i] = fileSplit[0]
		file2Split = re.split('\.', read2[i])
		sampleName2 = file2Split[0]
		assert sampleName[i] == sampleName2

		print("bwa mem {} {}/{} {}/{} > {}/sam/{}.sam".format(refSeq, read1Dir, read1[i], read2Dir, read2[i], out, sampleName[i]))
		os.system("bwa mem {} {}/{} {}/{} > {}/sam/{}.sam".format(refSeq, read1Dir, read1[i], read2Dir, read2[i], out, sampleName[i]))
		os.system("samtools view -b {}/sam/{}.sam > {}/bam/{}.bam".format(out, sampleName[i], out, sampleName[i]))
		os.system("rm {}/sam/{}.sam".format(out, sampleName[i]))
		os.system("samtools sort {}/bam/{}.bam -o {}/sort_bam/{}_sort.bam".format(out, sampleName[i], out, sampleName[i]))
		# os.system("rm {}/bam/{}.bam".format(out, sampleName[i]))
		os.system("samtools index -b {}/sort_bam/{}_sort.bam".format(out, sampleName[i]))

	return(sampleName)

def varCall(refSeq, samples, bq, out):
	anno = "FORMAT/AD,FORMAT/ADF,FORMAT/ADR,FORMAT/DP,FORMAT/SP,INFO/AD,INFO/ADF,INFO/ADR"
	size = len(samples)

	print("\nCalling variants\n")
	for i in range(size):
		os.system("bcftools mpileup -Ou --annotate {} -Q {} -d 300000 -f {} {}/sort_bam/{}_sort.bam \
		| bcftools call -Ov -m --skip-variants indels > \
		{}/vcf/{}.vcf".format(anno, bq, refSeq, out, samples[i], out, samples[i]))


#Main-------------------------------------------------------------------------------

numRef = 1
refs = snakemake.params["refs"]
pair1 = snakemake.params["trimmed"] + "/1"
pair2 = snakemake.params["trimmed"] + "/2"
baseQual = 20
out = snakemake.params["variant_out"]

reads1 = getReads(pair1)
reads2 = getReads(pair2)

# create the output directory if needed
if not os.path.exists(out):
	os.makedirs(out)
print("\nMinimum base quality set to {}\n".format(baseQual))
checkVersion()
makeOutDirs(out)
ref_fasta = getRef(refs)
print(ref_fasta)
samples = align(ref_fasta, reads1, reads2, pair1, pair2, out)
varCall(ref_fasta, samples, baseQual, out)