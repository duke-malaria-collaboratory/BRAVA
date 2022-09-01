configfile: "config.yaml"

# modifiable parameters
TARGET = config['target']
RECREATEREFFOLDER = config['recreateRefFolder']
CUTOFF = config['cutoff']
SEED = config['seed']
READ_DEPTH = config['read_depth']
PROPORTION = config['proportion']
HAPLOTYPE_LENGTH = config['haplotype_length']
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
 		expand("{out}/fastq/{target}/all_samples/", out=OUT, target=TARGET),
 		expand("{out}/haplotype_output/{target}_trimAndFilterTable", out=OUT, target=TARGET),
		expand("{out}/haplotype_output/{target}_haplotype_table_censored_final_version.csv", out=OUT, target=TARGET),

rule clean_sequencing_reads:
	input:
		refs=expand("refs/{target}.fasta", target=TARGET),
		pair1=PAIR1,
		pair2=PAIR2,
		forward=expand("{forward}.fasta", forward=FOWARD),
		rev=expand("{rev}.fasta", rev=REV)
	output:
		#directory(expand("{out}/fastq/{target}/1/", out=OUT, target=TARGET)),
		#directory(expand("{out}/fastq/{target}/2/", out=OUT, target=TARGET)),
		out=directory(expand("{out}/fastq/{target}/all_samples/", out=OUT, target=TARGET))
	params:
		recreate_ref_folder=RECREATEREFFOLDER,
		folder=OUT,
		output=expand("{out}/fastq/{target}/all_samples", out=OUT, target=TARGET),
		forward_samples=expand("{out}/fastq/{target}/1", out=OUT, target=TARGET),
		reverse_samples=expand("{out}/fastq/{target}/2", out=OUT, target=TARGET),
		haplotype_output=expand("{out}/haplotype_output", out=OUT),
		out=expand("{out}/fastq/{target}/all_samples", out=OUT, target=TARGET)
	priority:
			10
	run:
		if params.recreate_ref_folder == True:
			shell("rm -r ref")
		shell("perl scripts/step1_splitSyncReadsMultiRef.pl 1 {input.refs} {params.folder} {input.pair1} {input.pair2} {input.forward} {input.rev}")
		shell("mkdir {params.output}")
		shell("mv {params.forward_samples}/*.fastq.gz {params.out} && mv {params.reverse_samples}/*.fastq.gz {params.out}")
		shell("rm -r {params.forward_samples} && rm -r {params.reverse_samples}")

rule call_haplotypes:
	input:
		expand("{out}/fastq/{target}/all_samples/", out=OUT, target=TARGET),
	output:
		trim_filter_table=expand("{out}/haplotype_output/{target}_trimAndFilterTable", out=OUT, target=TARGET),
		results=expand("{out}/haplotype_output/{target}_haplotypes.rds", out=OUT, target=TARGET),
		reads_table=expand("{out}/haplotype_output/{target}_trackReadsThroughPipeline.csv", out=OUT, target=TARGET)
	params:
		all_samples=expand("{out}/fastq/{target}/all_samples", out=OUT, target=TARGET),
		rscript="scripts/step2_dada2.R",
		cutoff = CUTOFF,
		seed = SEED,
	priority:
		30
	script:
		"{params.rscript}"

rule censor_haplotypes:
	input:
		expand("{out}/haplotype_output/{target}_haplotypes.rds", out=OUT, target=TARGET),
	output:
		precensored_haplotype_table=expand("{out}/haplotype_output/{target}_haplotype_table_precensored.csv", out=OUT, target=TARGET),
		snps_between_haps=expand("{out}/haplotype_output/{target}_snps_between_haps_within_samples.fasta", out=OUT, target=TARGET),
		unique_seqs=expand("{out}/haplotype_output/{target}_uniqueSeqs.fasta", out=OUT, target=TARGET),
		aligned_seqs=expand("{out}/haplotype_output/{target}_aligned_seqs.fasta", out=OUT, target=TARGET),
		final_censored=expand("{out}/haplotype_output/{target}_uniqueSeqs_final_censored.fasta", out=OUT, target=TARGET),
		final_haplotype_table=expand("{out}/haplotype_output/{target}_haplotype_table_censored_final_version.csv", out=OUT, target=TARGET),
	params:
		rscript="scripts/step3_haplotype_censoring.R",
		haplotypes=expand("{out}/haplotype_output/{target}_haplotypes.rds", out=OUT, target=TARGET),
		depth=READ_DEPTH,
		proportion=PROPORTION,
		length=HAPLOTYPE_LENGTH,
		ratio=READ_DEPTH_RATIO,
	priority:
		40
	script:
		"{params.rscript}"