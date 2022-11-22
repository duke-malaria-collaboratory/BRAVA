import pandas as pd

configfile: "config.yaml"

# modifiable parameters - change in config.yaml file
TARGET_TABLE = pd.read_table(config['target_file'])
TARGET = list(TARGET_TABLE.target.unique())
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
		# expand("{out}/fastqc_in", out=OUT, target=TARGET),
		# expand("{out}/trim", out=OUT, target=TARGET),
		# expand("{target}/{out}/fastq/all_samples", out=OUT, target=TARGET),
 		# expand("{target}/{out}/haplotype_output/{target}_{q_values}_trimAndFilterTable", out=OUT, target=TARGET, q_values=TRUNCQ_VALUES),
		# expand("{target}/{out}/haplotype_output/{target}_finalTrimAndFilterTable", out=OUT, target=TARGET),
		# expand("{target}/{out}/haplotype_output/{target}_trackReadsThroughPipeline.csv", out=OUT, target=TARGET),
		expand("{target}/{out}/haplotype_output/{target}_haplotype_table_censored_final_version.csv", out=OUT, target=TARGET),

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
		pyscript="scripts/step1_callFastqc.py",
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
		pyscript="scripts/step2_callTrimmomatic.py"
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
		pyscript="scripts/step3_synchronizeReads.py",
	script:
		"{params.pyscript}"

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
		rscript="scripts/step4_trim_and_filter.R",
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
		rscript="scripts/step5_optimize_reads.R",
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
		rscript="scripts/step6_call_haplotypes.R",
	script:
		"{params.rscript}"

rule censor_haplotypes:
	input:
		"{target}/{out}/haplotype_output/{target}_haplotypes.rds",
	output:
		precensored_haplotype_table="{target}/{out}/haplotype_output/{target}_haplotype_table_precensored.csv",
		snps_between_haps="{target}/{out}/haplotype_output/{target}_snps_between_haps_within_samples.fasta",
		unique_seqs="{target}/{out}/haplotype_output/{target}_uniqueSeqs.fasta",
		aligned_seqs="{target}/{out}/haplotype_output/{target}_aligned_seqs.fasta",
		final_censored="{target}/{out}/haplotype_output/{target}_uniqueSeqs_final_censored.fasta",
		final_haplotype_table="{target}/{out}/haplotype_output/{target}_haplotype_table_censored_final_version.csv",
	params:
		haplotypes="{target}/{out}/haplotype_output/{target}_haplotypes.rds",
		depth=READ_DEPTH,
		proportion=PROPORTION,
		length=lambda wildcards: list(TARGET_TABLE.length[TARGET_TABLE.target == wildcards.target]),
		ratio=READ_DEPTH_RATIO,
		rscript="scripts/step7_haplotype_censoring.R",
	script:
		"{params.rscript}"