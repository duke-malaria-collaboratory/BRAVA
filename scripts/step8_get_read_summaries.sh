#!/usr/bin/env bash

echo "getting trim summaries"

ls ${snakemake_params[pre_trim_summaries]} > ${snakemake_output[trim_summary_names]}
cat ${snakemake_params[all_pre_trim_summaries]} > ${snakemake_output[trim_summaries]}

echo "getting pre-filt"

for i in $(ls ${snakemake_params[all_fastq_files]}); do paste <(echo $i) <(zcat $i | grep -c "^+$"); done > ${snakemake_output[pre_filt_fastq_counts]}

echo "getting filt"

for i in $(ls ${snakemake_params[all_filtered_files]}); do paste <(echo $i) <(zcat $i | grep -c "^+$"); done > ${snakemake_output[filt_fastq_counts]}
