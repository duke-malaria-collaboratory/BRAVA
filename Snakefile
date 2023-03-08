configfile: "config.yaml"

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
HAPLOTYPE_PAIR1 = config['haplotype_pair1']
HAPLOTYPE_PAIR2 = config['haplotype_pair2']
HAPLOTYPE_OUT = config['haplotype_out']
VARIANT_PAIR1 = config['variant_pair1']
VARIANT_PAIR2 = config['variant_pair2']
VARIANT_OUT = config['variant_out']

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
		# expand("{out}/multiqc_report.html", out=OUT),
		#expand("{out}/long_summary.csv", out=OUT),
		expand("{out}/dr_depths_freqs.csv", out=VARIANT_OUT) if (len(VARIANT_CALLING) > 0) else [],
		expand("{out}/long_summary.csv", out=HAPLOTYPE_OUT) if (len(HAPLOTYPE_CALLING) > 0) else [],

if (len(VARIANT_CALLING) > 0):
	rule call_fastqc_variant:
		input:
			pair1=VARIANT_PAIR1,
			pair2=VARIANT_PAIR2,
		output:
			variant_out=directory(expand("{variant_out}/fastqc_in/", variant_out=VARIANT_OUT)),
		params:
			out=VARIANT_OUT,
			pair1=VARIANT_PAIR1,
			pair2=VARIANT_PAIR2,
			pyscript="scripts/step1_call_fastqc.py",
		script:
			"{params.pyscript}"

	rule call_trimmomatic_variant:
		input:
			expand("{variant_out}/fastqc_in/", variant_out=VARIANT_OUT),
		output:
			out=directory(expand("{variant_out}/trim/", variant_out=VARIANT_OUT)),
		params:
			out=VARIANT_OUT,
			pair1=VARIANT_PAIR1,
			pair2=VARIANT_PAIR2,
			forward=expand("{forward}.fasta", forward=FOWARD),
			rev=expand("{rev}.fasta", rev=REV),
			pyscript="scripts/step2_call_trimmomatic.py"
		script:
			"{params.pyscript}"

	rule variant_calling: # calling variants
		input:
			"{variant_out}/trim/",
		output:
			variant_out=directory("{target}/{variant_out}/vcf"),
		params:
			variant_out="{target}/{variant_out}",
			refs="refs/{target}",
			trimmed="{variant_out}/trim",
			target="{target}",
			pyscript="scripts/step4a_variant_calling.py",
		script:
			"{params.pyscript}"
		
	rule analyze_vcf:
		input:
			"{target}/{variant_out}/vcf",
		output:
			"{target}/{variant_out}/DR_mutations.csv",
		params:
			table=VARIANT_TABLE,
			vcf="{target}/{variant_out}/vcf",
			target="{target}",
			pyscript="scripts/step4b_analyze_VCF.py",
		script:
			"{params.pyscript}"
		
	rule combine_vcf_analysis:
		input:
			expand("{target}/{variant_out}/DR_mutations.csv", target=VARIANT_CALLING, variant_out=VARIANT_OUT),
		output:
			dr_depths="{variant_out}/dr_depths_freqs.csv",
			dr_freqs="{variant_out}/dr_freqs.csv"
		params:
			targets=VARIANT_CALLING,
			variant_out=VARIANT_OUT,
			table=VARIANT_TABLE,
			rscript="scripts/step4c_combine_analysis.R",
		script:
			"{params.rscript}"

if (len(HAPLOTYPE_CALLING) > 0):
	rule call_fastqc_haplotype:
		input:
			pair1=HAPLOTYPE_PAIR1,
			pair2=HAPLOTYPE_PAIR2,
		output:
			out=directory(expand("{haplotype_out}/fastqc_in/", haplotype_out=HAPLOTYPE_OUT)),
		params:
			out=HAPLOTYPE_OUT,
			pair1=HAPLOTYPE_PAIR1,
			pair2=HAPLOTYPE_PAIR2,
			pyscript="scripts/step1_call_fastqc.py",
		script:
			"{params.pyscript}"

	rule call_trimmomatic_haplotype:
		input:
			expand("{haplotype_out}/fastqc_in/", haplotype_out=HAPLOTYPE_OUT)
		output:
			out=directory(expand("{haplotype_out}/trim/", haplotype_out=HAPLOTYPE_OUT)),
		params:
			out=HAPLOTYPE_OUT,
			pair1=HAPLOTYPE_PAIR1,
			pair2=HAPLOTYPE_PAIR2,
			forward=expand("{forward}.fasta", forward=FOWARD),
			rev=expand("{rev}.fasta", rev=REV),
			pyscript="scripts/step2_call_trimmomatic.py"
		script:
			"{params.pyscript}"

	rule generate_multiqc_report:
		input:
			"{haplotype_out}/trim",
		output:
			"{haplotype_out}/multiqc_report.html",
		params:
			haplotype_out="{haplotype_out}",
		shell:
			"multiqc . --outdir {params.haplotype_out}"
	
	rule get_marker_lengths:
		input:
			"{haplotype_out}/multiqc_report.html"
		output:
			"{haplotype_out}/marker_lengths.csv",
		params:
			refs=REFS,
			primers=expand("{forward}.fasta", forward=FOWARD),
			pyscript="scripts/step3_get_marker_lengths.py",
		script:
			"{params.pyscript}"

	rule synchronize_reads: # calling haplotypes
		input:
			"{haplotype_out}/trim",
		output:
			out=directory("{target}/{haplotype_out}/fastq/all_samples/"),
		params:
			out="{target}/{haplotype_out}",
			refs="refs/{target}/{target}.fasta",
			all_samples="{target}/{haplotype_out}/fastq/all_samples",
			forward_samples="{target}/{haplotype_out}/fastq/1",
			reverse_samples="{target}/{haplotype_out}/fastq/2",
			trimmed="{haplotype_out}/trim",
			target="{target}",
			pyscript="scripts/step4_synchronize_reads.py",
		script:
			"{params.pyscript}"

	rule trim_and_filter:
		input:
			"{target}/{haplotype_out}/fastq/all_samples/",
		output:
			trim_filter_table="{target}/{haplotype_out}/haplotype_output/{target}_{q_values}_trimAndFilterTable",
		params:
			all_samples="{target}/{haplotype_out}/fastq/all_samples",
			read_count="{target}/{haplotype_out}/haplotype_output/{target}_read_count",
			q_values="{q_values}",
			haplotype_output="{target}/{haplotype_out}/haplotype_output",
			rscript="scripts/step5_trim_and_filter.R",
		script:
			"{params.rscript}"

	rule optimize_reads:
		input:
			trim_filter_table=expand("{{target}}/{{haplotype_out}}/haplotype_output/{{target}}_{q_values}_trimAndFilterTable", q_values=TRUNCQ_VALUES),
			q_trim_filter=expand("{{target}}/{{haplotype_out}}/haplotype_output/{{target}}_{q_values}_trimAndFilterTable", q_values=TRUNCQ_VALUES),
		output:
			max_read_count="{target}/{haplotype_out}/haplotype_output/{target}_max_read_count",
			final_trim_filter_table="{target}/{haplotype_out}/haplotype_output/{target}_finalTrimAndFilterTable",
			final_q_value="{target}/{haplotype_out}/haplotype_output/{target}_final_q_value",
		params:
			out="{haplotype_out}",
			target="{target}",
			trim_filter_path="{target}/{haplotype_out}/haplotype_output/{target}",
			read_count="{target}/{haplotype_out}/haplotype_output/{target}_read_count",
			max_read_count="{target}/{haplotype_out}/haplotype_output/{target}_max_read_count",
			all_samples="{target}/{haplotype_out}/fastq/all_samples",
			final_filtered="{target}/{haplotype_out}/fastq/all_samples/final_filtered",
			rscript="scripts/step6_optimize_reads.R",
		script:
			"{params.rscript}"

	rule call_haplotypes:
		input:
			"{target}/{haplotype_out}/haplotype_output/{target}_finalTrimAndFilterTable",
		output:
			results="{target}/{haplotype_out}/haplotype_output/{target}_haplotypes.rds",
			reads_table="{target}/{haplotype_out}/haplotype_output/{target}_trackReadsThroughPipeline.csv",
		params:
			all_samples="{target}/{haplotype_out}/fastq/all_samples",
			trim_filter_table="{target}/{haplotype_out}/haplotype_output/{target}_finalTrimAndFilterTable",
			cutoff=CUTOFF,
			seed=SEED,
			rscript="scripts/step7_call_haplotypes.R",
		script:
			"{params.rscript}"

	rule censor_haplotypes:
		input:
			input_file="{target}/{haplotype_out}/haplotype_output/{target}_haplotypes.rds",
			lengths="{haplotype_out}/marker_lengths.csv",
		output:
			precensored_haplotype_table="{target}/{haplotype_out}/haplotype_output/{target}_haplotype_table_precensored.csv",
			snps_between_haps="{target}/{haplotype_out}/haplotype_output/{target}_snps_between_haps_within_samples.fasta",
			unique_seqs="{target}/{haplotype_out}/haplotype_output/{target}_uniqueSeqs.fasta",
			aligned_seqs="{target}/{haplotype_out}/haplotype_output/{target}_aligned_seqs.fasta",
			final_censored="{target}/{haplotype_out}/haplotype_output/{target}_uniqueSeqs_final_censored.fasta",
			final_haplotype_table="{target}/{haplotype_out}/haplotype_output/{target}_haplotype_table_censored_final_version.csv",
		params:
			target="{target}",
			haplotypes="{target}/{haplotype_out}/haplotype_output/{target}_haplotypes.rds",
			depth=READ_DEPTH,
			proportion=PROPORTION,
			marker_lengths="{haplotype_out}/marker_lengths.csv",
			ratio=READ_DEPTH_RATIO,
			rscript="scripts/step8_haplotype_censoring.R",
		script:
			"{params.rscript}"

	rule get_read_summaries:
		input: 
			expand("{target}/{haplotype_out}/haplotype_output/{target}_uniqueSeqs_final_censored.fasta", target=HAPLOTYPE_CALLING, haplotype_out=HAPLOTYPE_OUT),
		output:
			trim_summary_names="{haplotype_out}/trim_summary_names.txt",
			trim_summaries="{haplotype_out}/trim_summaries.txt",
			pre_filt_fastq_counts="{haplotype_out}/pre-filt_fastq_read_counts.txt",
			filt_fastq_counts="{haplotype_out}/filt_fastq_read_counts.txt",
		params:
			pre_trim_summaries="{haplotype_out}/trim/summaries",
			all_pre_trim_summaries="{haplotype_out}/trim/summaries/*",
			all_fastq_files="*/{haplotype_out}/fastq/all_samples/*fastq.gz",
			all_filtered_files="*/{haplotype_out}/fastq/all_samples/final_filtered/*",
		script:
			"scripts/step9_get_read_summaries.sh"

	rule create_summaries:
		input:
			trim_summaries="{haplotype_out}/trim_summaries.txt",
		output:
			long_summary="{haplotype_out}/long_summary.csv",
			wide_summary="{haplotype_out}/wide_summary.csv",
		params:
			trim_summary_names="{haplotype_out}/trim_summary_names.txt",
			trim_summaries="{haplotype_out}/trim_summaries.txt",
			pre_filt_fastq_counts="{haplotype_out}/pre-filt_fastq_read_counts.txt",
			filt_fastq_counts="{haplotype_out}/filt_fastq_read_counts.txt",
			all_filtered_files="/{haplotype_out}/fastq/all_samples/final_filtered/.*",
			rscript="scripts/step10_create_summaries.R",
		script:
			"{params.rscript}"
