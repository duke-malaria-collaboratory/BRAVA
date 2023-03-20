configfile: "config.yaml"

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
TRIM_FILTER_OUT = config['trim_filter_out']

# paths
REFS = config['refs']
FOWARD = config['forward']
REV = config['rev']
VARIANT_TABLE = config['variant_table']

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
		expand("{out}/multiqc_report.html", out=TRIM_FILTER_OUT),
		#expand("{out}/long_summary.csv", out=OUT),
		expand("{out}/dr_depths_freqs.csv", out=TRIM_FILTER_OUT) if (len(VARIANT_CALLING) > 0) else [],
		expand("{out}/long_summary.csv", out=TRIM_FILTER_OUT) if (len(HAPLOTYPE_CALLING) > 0) else [],

rule call_fastqc:
	input:
		pair1=PAIR1,
		pair2=PAIR2,
	output:
		trim_filter_out=directory(expand("{trim_filter_out}/fastqc_in/", trim_filter_out=TRIM_FILTER_OUT)),
	params:
		out=TRIM_FILTER_OUT,
		pair1=PAIR1,
		pair2=PAIR2,
		pyscript="scripts/step1_call_fastqc.py",
	script:
		"{params.pyscript}"

rule call_trimmomatic:
	input:
		expand("{trim_filter_out}/fastqc_in/", trim_filter_out=TRIM_FILTER_OUT),
	output:
		out=directory(expand("{trim_filter_out}/trim/", trim_filter_out=TRIM_FILTER_OUT)),
	params:
		out=TRIM_FILTER_OUT,
		pair1=PAIR1,
		pair2=PAIR2,
		forward=expand("{forward}.fasta", forward=FOWARD),
		rev=expand("{rev}.fasta", rev=REV),
		pyscript="scripts/step2_call_trimmomatic.py"
	script:
		"{params.pyscript}"

rule generate_multiqc_report:
	input:
		"{trim_filter_out}/trim",
	output:
		"{trim_filter_out}/multiqc_report.html",
	params:
		trim_filter_out="{trim_filter_out}",
	shell:
		"multiqc . --outdir {params.trim_filter_out}"

# variant calling
if (len(VARIANT_CALLING) > 0):
	rule variant_calling:
		input:
			expand("{trim_filter_out}/trim/", trim_filter_out=TRIM_FILTER_OUT),
		output:
			variant_out=directory("{target}/out/vcf"),
		params:
			variant_out="{target}/out",
			refs="refs/{target}",
			trimmed=expand("{trim_filter_out}/trim", trim_filter_out=TRIM_FILTER_OUT),
			target="{target}",
			pyscript="scripts/step4a_variant_calling.py",
		script:
			"{params.pyscript}"
		
	rule analyze_vcf:
		input:
			"{target}/out/vcf",
		output:
			"{target}/out/DR_mutations.csv",
		params:
			table=VARIANT_TABLE,
			vcf="{target}/out/vcf",
			target="{target}",
			pyscript="scripts/step4b_analyze_VCF.py",
		script:
			"{params.pyscript}"
		
	rule combine_vcf_analysis:
		input:
			expand("{target}/out/DR_mutations.csv", target=VARIANT_CALLING),
		output:
			dr_depths="{trim_filter_out}/dr_depths_freqs.csv",
		params:
			targets=VARIANT_CALLING,
			trim_filter_out=TRIM_FILTER_OUT,
			table=VARIANT_TABLE,
			DR_mutations=expand("{target}/out/DR_mutations.csv", target=VARIANT_CALLING),
			rscript="scripts/step4c_combine_analysis.R",
		script:
			"{params.rscript}"

# haplotype calling
if (len(HAPLOTYPE_CALLING) > 0):
	rule get_marker_lengths:
		input:
			expand("{trim_filter_out}/trim", trim_filter_out=TRIM_FILTER_OUT),
		output:
			"{trim_filter_out}/marker_lengths.csv",
		params:
			refs=REFS,
			primers=expand("{forward}.fasta", forward=FOWARD),
			pyscript="scripts/step3_get_marker_lengths.py",
		script:
			"{params.pyscript}"

	rule synchronize_reads:
		input:
			expand("{trim_filter_out}/trim", trim_filter_out=TRIM_FILTER_OUT),
		output:
			out=directory("{target}/out/fastq/all_samples/"),
		params:
			out="{target}/out",
			refs="refs/{target}/{target}.fasta",
			all_samples="{target}/out/fastq/all_samples",
			forward_samples="{target}/out/fastq/1",
			reverse_samples="{target}/out/fastq/2",
			trimmed=expand("{trim_filter_out}/trim", trim_filter_out=TRIM_FILTER_OUT),
			target="{target}",
			pyscript="scripts/step4_synchronize_reads.py",
		script:
			"{params.pyscript}"

	rule trim_and_filter:
		input:
			"{target}/out/fastq/all_samples/",
		output:
			trim_filter_table="{target}/out/haplotype_output/{target}_{q_values}_trimAndFilterTable",
		params:
			all_samples="{target}/out/fastq/all_samples",
			read_count="{target}/out/haplotype_output/{target}_read_count",
			q_values="{q_values}",
			haplotype_output="{target}/out/haplotype_output",
			rscript="scripts/step5_trim_and_filter.R",
		script:
			"{params.rscript}"

	rule optimize_reads:
		input:
			trim_filter_table=expand("{{target}}/out/haplotype_output/{{target}}_{q_values}_trimAndFilterTable", q_values=TRUNCQ_VALUES),
			q_trim_filter=expand("{{target}}/out/haplotype_output/{{target}}_{q_values}_trimAndFilterTable", q_values=TRUNCQ_VALUES),
		output:
			max_read_count="{target}/out/haplotype_output/{target}_max_read_count",
			final_trim_filter_table="{target}/out/haplotype_output/{target}_finalTrimAndFilterTable",
			final_q_value="{target}/out/haplotype_output/{target}_final_q_value",
		params:
			out="out",
			target="{target}",
			trim_filter_path="{target}/out/haplotype_output/{target}",
			read_count="{target}/out/haplotype_output/{target}_read_count",
			max_read_count="{target}/out/haplotype_output/{target}_max_read_count",
			all_samples="{target}/out/fastq/all_samples",
			final_filtered="{target}/out/fastq/all_samples/final_filtered",
			rscript="scripts/step6_optimize_reads.R",
		script:
			"{params.rscript}"

	rule call_haplotypes:
		input:
			"{target}/out/haplotype_output/{target}_finalTrimAndFilterTable",
		output:
			results="{target}/out/haplotype_output/{target}_haplotypes.rds",
			reads_table="{target}/out/haplotype_output/{target}_trackReadsThroughPipeline.csv",
		params:
			all_samples="{target}/out/fastq/all_samples",
			trim_filter_table="{target}/out/haplotype_output/{target}_finalTrimAndFilterTable",
			cutoff=CUTOFF,
			seed=SEED,
			rscript="scripts/step7_call_haplotypes.R",
		script:
			"{params.rscript}"

	rule censor_haplotypes:
		input:
			input_file="{target}/out/haplotype_output/{target}_haplotypes.rds",
			lengths=expand("{trim_filter_out}/marker_lengths.csv", trim_filter_out=TRIM_FILTER_OUT),
		output:
			precensored_haplotype_table="{target}/out/haplotype_output/{target}_haplotype_table_precensored.csv",
			snps_between_haps="{target}/out/haplotype_output/{target}_snps_between_haps_within_samples.fasta",
			unique_seqs="{target}/out/haplotype_output/{target}_uniqueSeqs.fasta",
			aligned_seqs="{target}/out/haplotype_output/{target}_aligned_seqs.fasta",
			final_censored="{target}/out/haplotype_output/{target}_uniqueSeqs_final_censored.fasta",
			final_haplotype_table="{target}/out/haplotype_output/{target}_haplotype_table_censored_final_version.csv",
		params:
			target="{target}",
			haplotypes="{target}/out/haplotype_output/{target}_haplotypes.rds",
			depth=READ_DEPTH,
			proportion=PROPORTION,
			marker_lengths=expand("{trim_filter_out}/marker_lengths.csv", trim_filter_out=TRIM_FILTER_OUT),
			ratio=READ_DEPTH_RATIO,
			rscript="scripts/step8_haplotype_censoring.R",
		script:
			"{params.rscript}"

	rule get_read_summaries:
		input: 
			expand("{target}/out/haplotype_output/{target}_uniqueSeqs_final_censored.fasta", target=HAPLOTYPE_CALLING),
		output:
			trim_summary_names="{trim_filter_out}/trim_summary_names.txt",
			trim_summaries="{trim_filter_out}/trim_summaries.txt",
			pre_filt_fastq_counts="{trim_filter_out}/pre-filt_fastq_read_counts.txt",
			filt_fastq_counts="{trim_filter_out}/filt_fastq_read_counts.txt",
		params:
			pre_trim_summaries="{trim_filter_out}/trim/summaries",
			all_pre_trim_summaries="{trim_filter_out}/trim/summaries/*",
			all_fastq_files="*/out/fastq/all_samples/*fastq.gz",
			all_filtered_files="*/out/fastq/all_samples/final_filtered/*",
		script:
			"scripts/step9_get_read_summaries.sh"

# keep censored, precensored, remove q value stuff and .rds
	rule create_summaries:
		input:
			trim_summaries="{trim_filter_out}/trim_summaries.txt",
		output:
			long_summary="{trim_filter_out}/long_summary.csv",
			wide_summary="{trim_filter_out}/wide_summary.csv",
		params:
			trim_summary_names="{trim_filter_out}/trim_summary_names.txt",
			trim_summaries="{trim_filter_out}/trim_summaries.txt",
			pre_filt_fastq_counts="{trim_filter_out}/pre-filt_fastq_read_counts.txt",
			filt_fastq_counts="{trim_filter_out}/filt_fastq_read_counts.txt",
			all_filtered_files="/{trim_filter_out}/fastq/all_samples/final_filtered/.*",
			rscript="scripts/step10_create_summaries.R",
		script:
			"{params.rscript}"
