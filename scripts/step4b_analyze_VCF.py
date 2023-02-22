#!/usr/bin/env python

import pandas as pd
import os
import io
from os import listdir

"""
This script is called as a single sample submission or Slurm array job, and returns
a table of drug-resistance-associated mutations among amplicon vcf files.

*************************************************************************************

command: python analyzeVCF.py
"""

#Functions--------------------------------------------------------------------------
def get_samples(path):
    mypath = path
    sampleList = []
    sampleFileList = []
    for file in os.listdir(mypath):
        if file.endswith(".vcf"):
            sample = file.split('.')[0]
            sampleList.append(sample)
            sampleFileList.append(os.path.join(mypath, file))
    sampleFileList.sort()
    sampleList.sort()
    return (sampleFileList, sampleList)

def read_vcf(path):
    with open(path, 'r') as f:
        lines = [l for l in f if not l.startswith('##')]
    return pd.read_csv(
        io.StringIO(''.join(lines)),
        dtype={'#CHROM': str, 'POS': int, 'ID': str, 'REF': str, 'ALT': str,
               'QUAL': str, 'FILTER': str, 'INFO': str},
        sep='\t'
    ).rename(columns={'#CHROM': 'CHROM'})

def get_format(df, sample):
    #get list of info fields for each variant split on ';'
    IDs = df.iloc[:,8].str.split(':')
    form = df.iloc[:,9].str.split(':')
    varIDs = []
    #identify variant FORMAT field
    for i in range(len(IDs)):
        if (len(IDs[i]) == 7):
            varIDs = IDs[i]
            break
    # identify positions of DP and AD in variant format field
    for i in range(len(varIDs)):
        if varIDs[i] == 'DP':
            dp_loc = i
        elif varIDs[i] == 'AD':
            ad_loc = i

    #length of list is equal to number of variants in sample
    numVars = len(form)
    myDict = {}
    #iterate over variants called in sample
    for i in range(numVars):
        pos = df['POS'][i]
        ref = df['REF'][i]
        alt = df['ALT'][i]

        if len(form[i]) < 7:
            dp = form[i][1]
            ad = 0
            alt_freq = 0
        elif len(form[i]) == 7:
            dp = form[i][dp_loc]
            ad = form[i][ad_loc].split(',')[1]
            alt_freq = int(ad)/int(dp)

        myDict[pos] = {}
        myDict[pos]['REF'] = ref
        myDict[pos]['ALT'] = alt
        myDict[pos]['ALT_DEPTH'] = ad
        myDict[pos]['DEPTH'] = dp
        myDict[pos]['ALT_FREQ'] = alt_freq
    return myDict

def dict_to_pd(dic):
    df = pd.DataFrame.from_dict({(i,j,k): dic[i][j][k] 
        for i in dic.keys() 
            for j in dic[i].keys()
                for k in dic[i][j].keys()},
        orient='index')
    return df

def VCF_dataframe(targetName, files, names):

    target_dict = {}
    target = targetName
    target_dict[target] = {}
    for i in range(len(files)):
        sampleName = names[i]
        vcf = read_vcf(files[i])
        sample = get_format(vcf, names[i])
        target_dict[target][sampleName] = sample
    df = dict_to_pd(target_dict)
    df = df.reset_index()
    df = df.rename(columns={"level_0": "Target", "level_1": "Sample", "level_2":"POS"})
    return df

#Main-------------------------------------------------------------------------------
if __name__ == "__main__":
    
    description = "Organizes P. falciparum variants associated with drug resistance from amplicon deep sequencing" \

    refTable = pd.read_csv(snakemake.params["table"])
    refTable['ID'] = refTable['Target'] + refTable['Pos'].astype(str)
    refs = refTable[['ID']]

    files, samples = get_samples(snakemake.params["vcf"])
    df = VCF_dataframe(snakemake.params["target"], files, samples)
    allSNPs = df
    
    # dhfr_files, dhfr_samples = get_samples("/data/taylorlab/jws48/variantCall/dhfr/vcf/")
    # dhfr = VCF_dataframe("dhfr", dhfr_files, dhfr_samples)
    # allSNPs = dhfr
    
    # dhps_files, dhps_samples = get_samples("/data/taylorlab/jws48/variantCall/dhps/vcf/")
    # dhps = VCF_dataframe("dhps", dhps_files, dhps_samples)
    # allSNPs = allSNPs.append(dhps)

    # K13_files, K13_samples = get_samples("/data/taylorlab/jws48/variantCall/K13/vcf/")
    # K13 = VCF_dataframe("K13", K13_files, K13_samples)
    # allSNPs = allSNPs.append(K13)

    # pfmdr1_files, pfmdr1_samples = get_samples("/data/taylorlab/jws48/variantCall/pfmdr1/vcf/")
    # pfmdr1 = VCF_dataframe("pfmdr1", pfmdr1_files, pfmdr1_samples)
    # allSNPs = allSNPs.append(pfmdr1)

    # pfcrt_files, pfcrt_samples = get_samples("/data/taylorlab/jws48/variantCall/pfcrt/vcf/")
    # pfcrt = VCF_dataframe("pfcrt", pfcrt_files, pfcrt_samples)
    # allSNPs = allSNPs.append(pfcrt)
    allSNPs['ID'] = allSNPs['Target'] + allSNPs['POS'].astype(str)
    
    DR_mutations = allSNPs.merge(refs, how='right', left_on='ID', right_on='ID').drop(columns='ID')
    DR_mutations = DR_mutations.dropna().astype({"POS": int})
    DR_mutations.to_csv(snakemake.output[0])