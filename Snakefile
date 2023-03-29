configfile: "config/config.yaml"

# work on haplotype/variant calling stuff getting fixed first, update readme, then clean up output files

# modifiable parameters - change in config.yaml file
#TARGET = config['target']
VARIANT_CALLING = config['variant_calling_targets']
HAPLOTYPE_CALLING = config['haplotype_calling_targets']
TRUNCQ_VALUES = config['truncQ_values']
CUTOFF = config['cutoff']
SEED = config['seed']
READ_DEPTH = config['read_depth']
PROPORTION = config['proportion']
READ_DEPTH_RATIO = config['read_depth_ratio']

# input & output
PAIR1 = config['pair1']
PAIR2 = config['pair2']

# paths
REFS = config['refs']
FOWARD = config['forward']
REV = config['rev']
VARIANT_TABLE = config['variant_table']
ROOT = config['root']

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
		# expand("{out}/dr_depths_freqs.csv", out=OUT),
		#expand("{target}/out/fastq/all_samples/", target=HAPLOTYPE_CALLING),
		expand("{root}/output/multiqc_report.html", root=ROOT),
		expand("{root}/output/variant_output/dr_depths_freqs.csv", root=ROOT) if (len(VARIANT_CALLING) > 0) else [],
		expand("{root}/output/haplotype_output/long_summary.csv", root=ROOT) if (len(HAPLOTYPE_CALLING) > 0) else [],

rule call_fastqc:
	input:
		pair1=PAIR1,
		pair2=PAIR2,
	output:
		fastq_out=directory("{root}/output/fastqc_out/"),
	params:
		out="{root}/output",
		pair1=PAIR1,
		pair2=PAIR2,
		pyscript="scripts/step1_call_fastqc.py",
	script:
		"{params.pyscript}"

rule call_trimmomatic:
	input:
		"{root}/output/fastqc_out/",
	output:
		trim_filter_out=directory("{root}/output/trimmed_reads/"),
	params:
		out="{root}/output/trimmed_reads",
		pair1=PAIR1,
		pair2=PAIR2,
		forward=expand("{forward}.fasta", forward=FOWARD),
		rev=expand("{rev}.fasta", rev=REV),
		pyscript="scripts/step2_call_trimmomatic.py"
	script:
		"{params.pyscript}"

rule generate_multiqc_report:
	input:
		"{root}/output/trimmed_reads/",
	output:
		"{root}/output/multiqc_report.html",
	params:
		out="{root}/output",
	shell:
		"multiqc . --outdir {params.out}"

# variant calling
if (len(VARIANT_CALLING) > 0):
	rule variant_calling:
		input:
			"{root}/output/trimmed_reads/",
		output:
			variant_out=directory("{root}/output/variant_output/{target}/vcf"),
		params:
			variant_out="{root}/output/variant_output/{target}",
			refs=REFS + "/{target}",
			trimmed="{root}/output/trimmed_reads/",
			target="{target}",
			pyscript="scripts/step3_variant_calling.py",
		script:
			"{params.pyscript}"
		
	rule analyze_vcf:
		input:
			"{root}/output/variant_output/{target}/vcf",
		output:
			"{root}/output/variant_output/{target}/DR_mutations.csv",
		params:
			table=VARIANT_TABLE,
			vcf="{root}/output/variant_output/{target}/vcf",
			target="{target}",
			pyscript="scripts/step4_analyze_VCF.py",
		script:
			"{params.pyscript}"
		
	rule combine_vcf_analysis:
		input:
			expand("{{root}}/output/variant_output/{target}/DR_mutations.csv", target=VARIANT_CALLING),
		output:
			dr_depths="{root}/output/variant_output/dr_depths_freqs.csv",
		params:
			path="{root}/output/variant_output/",
			targets=VARIANT_CALLING,
			table=VARIANT_TABLE,
			DR_mutations=expand("{root}/output/variant_output/{target}/DR_mutations.csv", root=ROOT, target=VARIANT_CALLING),
			rscript="scripts/step5_combine_analysis.R",
		script:
			"{params.rscript}"

# haplotype calling
if (len(HAPLOTYPE_CALLING) > 0):
	rule get_marker_lengths:
		input:
			"{root}/output/trimmed_reads/",
		output:
			"{root}/output/marker_lengths.csv",
		params:
			refs=REFS,
			primers=expand("{forward}.fasta", forward=FOWARD),
			pyscript="scripts/step6_get_marker_lengths.py",
		script:
			"{params.pyscript}"

	rule synchronize_reads:
		input:
			"{root}/output/trimmed_reads/",
		output:
			out=directory("{root}/output/haplotype_output/{target}/bbsplit_out/"),
		params:
			out="{root}/output/haplotype_output/{target}/bbsplit_out",
			refs=REFS + "/{target}/{target}.fasta",
			mapped_reads="{root}/output/haplotype_output/{target}/bbsplit_out/mapped_reads",
			forward_samples="{root}/output/haplotype_output/{target}/bbsplit_out/1",
			reverse_samples="{root}/output/haplotype_output/{target}/bbsplit_out/2",
			trimmed="{root}/output/trimmed_reads",
			target="{target}",
			pyscript="scripts/step7_synchronize_reads.py",
		script:
			"{params.pyscript}"

	rule trim_and_filter:
		input:
			"{root}/output/haplotype_output/{target}/bbsplit_out/",
		output:
			trim_filter_table="{root}/output/haplotype_output/{target}/trim_filter_out/{target}_{q_values}_trim_and_filter_table",
		params:
			mapped_reads="{root}/output/haplotype_output/{target}/bbsplit_out/mapped_reads",
			read_count="{root}/output/haplotype_output/{target}/trim_filter_out/{target}_read_count",
			q_values="{q_values}",
			trim_filter_out="{root}/output/haplotype_output/{target}/trim_filter_out",
			rscript="scripts/step8_trim_and_filter.R",
		script:
			"{params.rscript}"

	rule optimize_reads:
		input:
			trim_filter_table=expand("{{root}}/output/haplotype_output/{{target}}/trim_filter_out/{{target}}_{q_values}_trim_and_filter_table", q_values=TRUNCQ_VALUES),
		output:
			max_read_count="{root}/output/haplotype_output/{target}/optimize_reads_out/{target}_max_read_count",
			final_trim_filter_table="{root}/output/haplotype_output/{target}/optimize_reads_out/{target}_final_trim_and_filter_table",
			final_q_value="{root}/output/haplotype_output/{target}/optimize_reads_out/{target}_final_q_value",
		params:
			out="out",
			target="{target}",
			trim_filter_path="{root}/output/haplotype_output/{target}/trim_filter_out/{target}",
			read_count="{root}/output/haplotype_output/{target}/trim_filter_out/{target}_read_count",
			max_read_count="{root}/output/haplotype_output/{target}/optimize_reads_out/{target}_max_read_count",
			mapped_reads="{root}/output/haplotype_output/{target}/bbsplit_out/mapped_reads",
			final_filtered="{root}/output/haplotype_output/{target}/bbsplit_out/mapped_reads/final_filtered",
			rscript="scripts/step9_optimize_reads.R",
		script:
			"{params.rscript}"

	rule call_haplotypes:
		input:
			"{root}/output/haplotype_output/{target}/optimize_reads_out/{target}_final_trim_and_filter_table",
		output:
			results="{root}/output/haplotype_output/{target}/{target}_haplotypes.rds",
			reads_table="{root}/output/haplotype_output/{target}/trim_filter_out/{target}_track_reads_through_pipeline.csv",
		params:
			mapped_reads="{root}/output/haplotype_output/{target}/bbsplit_out/mapped_reads",
			trim_filter_table="{root}/output/haplotype_output/{target}/optimize_reads_out/{target}_final_trim_and_filter_table",
			cutoff=CUTOFF,
			seed=SEED,
			rscript="scripts/step10_call_haplotypes.R",
		script:
			"{params.rscript}"

	rule censor_haplotypes:
		input:
			input_file="{root}/output/haplotype_output/{target}/{target}_haplotypes.rds",
			lengths="{root}/output/marker_lengths.csv",
		output:
			precensored_haplotype_table="{root}/output/haplotype_output/{target}/haplotypes/{target}_haplotype_table_precensored.csv",
			snps_between_haps="{root}/output/haplotype_output/{target}/sequences/{target}_snps_between_haps_within_samples.fasta",
			unique_seqs="{root}/output/haplotype_output/{target}/sequences/{target}_unique_seqs.fasta",
			aligned_seqs="{root}/output/haplotype_output/{target}/sequences/{target}_aligned_seqs.fasta",
			final_censored="{root}/output/haplotype_output/{target}/sequences/{target}_unique_seqs_final_censored.fasta",
			final_haplotype_table="{root}/output/haplotype_output/{target}/haplotypes/{target}_haplotype_table_censored_final_version.csv",
		params:
			target="{target}",
			haplotypes="{root}/output/haplotype_output/{target}/{target}_haplotypes.rds",
			depth=READ_DEPTH,
			proportion=PROPORTION,
			marker_lengths="{root}/output/marker_lengths.csv",
			ratio=READ_DEPTH_RATIO,
			rscript="scripts/step11_haplotype_censoring.R",
		script:
			"{params.rscript}"

	rule get_read_summaries:
		input: 
			expand("{{root}}/output/haplotype_output/{target}/sequences/{target}_unique_seqs_final_censored.fasta", target=HAPLOTYPE_CALLING),
		output:
			pre_filt_fastq_counts="{root}/output/haplotype_output/filtered_dada2/pre-filt_fastq_read_counts.txt",
			filt_fastq_counts="{root}/output/haplotype_output/filtered_dada2/filt_fastq_read_counts.txt",
		params:
			trim_summary_names="{root}/output/trimmed_reads/trim_summary_names.txt",
			trim_summaries="{root}/output/trimmed_reads/trim_summaries.txt",
			pre_trim_summaries="{root}/output/trimmed_reads/summaries",
			all_pre_trim_summaries="{root}/output/trimmed_reads/summaries/*",
			all_fastq_files="{root}/output/haplotype_output/*/bbsplit_out/mapped_reads/*fastq.gz",
			all_filtered_files="{root}/output/haplotype_output/*/bbsplit_out/mapped_reads/final_filtered/*",
		script:
			"scripts/step12_get_read_summaries.sh"

	rule create_summaries:
		input:
			pre_filt_fastq_counts="{root}/output/haplotype_output/filtered_dada2/pre-filt_fastq_read_counts.txt",
		output:
			long_summary="{root}/output/haplotype_output/long_summary.csv",
			wide_summary="{root}/output/haplotype_output/wide_summary.csv",
		params:
			trim_summary_names="{root}/output/trimmed_reads/trim_summary_names.txt",
			trim_summaries="{root}/output/trimmed_reads/trim_summaries.txt",
			pre_filt_fastq_counts="{root}/output/haplotype_output/filtered_dada2/pre-filt_fastq_read_counts.txt",
			filt_fastq_counts="{root}/output/haplotype_output/filtered_dada2/filt_fastq_read_counts.txt",
			all_filtered_files="{root}/output/haplotype_output//bbsplit_out/mapped_reads/final_filtered/.*",
			rscript="scripts/step13_create_summaries.R",
		script:
			"{params.rscript}"
