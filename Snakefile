configfile: "config.yaml"

# modifiable parameters - change in config.yaml file
TARGET = config['target']
TRUNCQ_VALUES = config['truncQ_values']
CUTOFF = config['cutoff']
SEED = config['seed']
READ_DEPTH = config['read_depth']
PROPORTION = config['proportion']
READ_DEPTH_RATIO = config['read_depth_ratio']

# paths
REFS = config['refs']
PAIR1 = config['pair1']
PAIR2 = config['pair2']
FOWARD = config['forward']
REV = config['rev']
OUT = config['out']

rule all:
	input:
		# expand("out/fastqc_in", out=OUT, target=TARGET),
		# expand("out/trim", out=OUT, target=TARGET),
		# expand("{target}/{out}/fastq/all_samples", out=OUT, target=TARGET),
 		# expand("{target}/{out}/haplotype_output/{target}_{q_values}_trimAndFilterTable", out=OUT, target=TARGET, q_values=TRUNCQ_VALUES),
		# expand("{target}/{out}/haplotype_output/{target}_finalTrimAndFilterTable", out=OUT, target=TARGET),
		# expand("{target}/{out}/haplotype_output/{target}_trackReadsThroughPipeline.csv", out=OUT, target=TARGET),
		# expand("out/trim_summary", out=OUT),
		# expand("{target}/{out}/haplotype_output/{target}_haplotype_table_censored_final_version.csv", out=OUT, target=TARGET),
		# expand("out/trim_summaries.txt", out=OUT),
		expand("{out}/multiqc_report.html", out=OUT),
		expand("{out}/long_summary.csv", out=OUT),

rule call_fastqc:
	input:
		pair1=PAIR1,
		pair2=PAIR2,
	output:
		out=directory("{out}/fastqc_in/"),
	params:
		out="{out}",
		pair1=PAIR1,
		pair2=PAIR2,
		pyscript="scripts/step1_call_fastqc.py",
	script:
		"{params.pyscript}"

rule call_trimmomatic:
	input:
		"{out}/fastqc_in/",
	output:
		out=directory("{out}/trim/"),
	params:
		out="{out}",
		pair1=PAIR1,
		pair2=PAIR2,
		forward=expand("{forward}.fasta", forward=FOWARD),
		rev=expand("{rev}.fasta", rev=REV),
		pyscript="scripts/step2_call_trimmomatic.py"
	script:
		"{params.pyscript}"

rule generate_multiqc_report:
	input:
		"{out}/trim",
	output:
		"{out}/multiqc_report.html",
	shell:
		"multiqc . --outdir out"
	
rule get_marker_lengths:
	input:
		"{out}/multiqc_report.html"
	output:
		"{out}/marker_lengths.csv",
	params:
		refs="refs",
		primers=expand("{forward}.fasta", forward=FOWARD),
		pyscript="scripts/step3_get_marker_lengths.py",
	script:
		"{params.pyscript}"
		
rule synchronize_reads:
	input:
		"{out}/trim/",
	output:
		out=directory("{target}/{out}/fastq/all_samples/"),
	params:
		out="{target}/{out}",
		refs="refs/{target}.fasta",
		all_samples="{target}/{out}/fastq/all_samples",
		forward_samples="{target}/{out}/fastq/1",
		reverse_samples="{target}/{out}/fastq/2",
		trimmed="{out}/trim",
		target="{target}",
		pyscript="scripts/step4_synchronizeReads.py",
	script:
		"{params.pyscript}"

# variant calling!
# then variantArray, then analyze VCF
# download the repo and try running it to see what the output files are

rule trim_and_filter:
	input:
		"{target}/{out}/fastq/all_samples/",
	output:
		trim_filter_table="{target}/{out}/haplotype_output/{target}_{q_values}_trimAndFilterTable",
	params:
		all_samples="{target}/{out}/fastq/all_samples",
		read_count="{target}/{out}/haplotype_output/{target}_read_count",
		q_values="{q_values}",
		haplotype_output="{target}/{out}/haplotype_output",
		rscript="scripts/step5_trim_and_filter.R",
	script:
		"{params.rscript}"

rule optimize_reads:
	input:
		trim_filter_table=expand("{target}/{out}/haplotype_output/{target}_{q_values}_trimAndFilterTable", out=OUT, target=TARGET, q_values=TRUNCQ_VALUES),
		q_trim_filter=expand("{target}/{out}/haplotype_output/{target}_{q_values}_trimAndFilterTable", out=OUT, target=TARGET, q_values=TRUNCQ_VALUES),
	output:
		max_read_count="{target}/{out}/haplotype_output/{target}_max_read_count",
		final_trim_filter_table="{target}/{out}/haplotype_output/{target}_finalTrimAndFilterTable",
		final_q_value="{target}/{out}/haplotype_output/{target}_final_q_value",
	params:
		out="{out}",
		target="{target}",
		read_count="{target}/{out}/haplotype_output/{target}_read_count",
		max_read_count="{target}/{out}/haplotype_output/{target}_max_read_count",
		all_samples="{target}/{out}/fastq/all_samples",
		final_filtered="{target}/{out}/fastq/all_samples/final_filtered",
		rscript="scripts/step6_optimize_reads.R",
	script:
		"{params.rscript}"

rule call_haplotypes:
	input:
		"{target}/{out}/haplotype_output/{target}_finalTrimAndFilterTable",
	output:
		results="{target}/{out}/haplotype_output/{target}_haplotypes.rds",
		reads_table="{target}/{out}/haplotype_output/{target}_trackReadsThroughPipeline.csv",
	params:
		all_samples="{target}/{out}/fastq/all_samples",
		trim_filter_table="{target}/{out}/haplotype_output/{target}_finalTrimAndFilterTable",
		cutoff=CUTOFF,
		seed=SEED,
		rscript="scripts/step7_call_haplotypes.R",
	script:
		"{params.rscript}"

rule censor_haplotypes:
	input:
		input_file="{target}/{out}/haplotype_output/{target}_haplotypes.rds",
		lengths="{out}/marker_lengths.csv",
	output:
		precensored_haplotype_table="{target}/{out}/haplotype_output/{target}_haplotype_table_precensored.csv",
		snps_between_haps="{target}/{out}/haplotype_output/{target}_snps_between_haps_within_samples.fasta",
		unique_seqs="{target}/{out}/haplotype_output/{target}_uniqueSeqs.fasta",
		aligned_seqs="{target}/{out}/haplotype_output/{target}_aligned_seqs.fasta",
		final_censored="{target}/{out}/haplotype_output/{target}_uniqueSeqs_final_censored.fasta",
		final_haplotype_table="{target}/{out}/haplotype_output/{target}_haplotype_table_censored_final_version.csv",
	params:
		target="{target}",
		haplotypes="{target}/{out}/haplotype_output/{target}_haplotypes.rds",
		depth=READ_DEPTH,
		proportion=PROPORTION,
		marker_lengths="{out}/marker_lengths.csv",
		ratio=READ_DEPTH_RATIO,
		rscript="scripts/step8_haplotype_censoring.R",
	script:
		"{params.rscript}"

rule get_read_summaries:
	input: 
		expand("{target}/{out}/haplotype_output/{target}_uniqueSeqs_final_censored.fasta", target=TARGET, out=OUT),
	output:
		trim_summary_names="{out}/trim_summary_names.txt",
		trim_summaries="{out}/trim_summaries.txt",
		pre_filt_fastq_counts="{out}/pre-filt_fastq_read_counts.txt",
		filt_fastq_counts="{out}/filt_fastq_read_counts.txt",
	params:
		pre_trim_summaries="{out}/trim/summaries",
		all_pre_trim_summaries="{out}/trim/summaries/*",
		all_fastq_files="*/{out}/fastq/all_samples/*fastq.gz",
		all_filtered_files="*/{out}/fastq/all_samples/final_filtered/*",
	script:
		"scripts/step9_get_read_summaries.sh"

rule create_summaries:
	input:
		trim_summaries="{out}/trim_summaries.txt",
	output:
		long_summary="{out}/long_summary.csv",
		wide_summary="{out}/wide_summary.csv",
	params:
		trim_summary_names="{out}/trim_summary_names.txt",
		trim_summaries="{out}/trim_summaries.txt",
		pre_filt_fastq_counts="{out}/pre-filt_fastq_read_counts.txt",
		filt_fastq_counts="{out}/filt_fastq_read_counts.txt",
		all_filtered_files="/{out}/fastq/all_samples/final_filtered/.*",
		rscript="scripts/step10_create_summaries.R",
	script:
		"{params.rscript}"
