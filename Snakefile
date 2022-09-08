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
		refs="refs/{target}.fasta",
		pair1=PAIR1,
		pair2=PAIR2,
		forward=expand("{forward}.fasta", forward=FOWARD),
		rev=expand("{rev}.fasta", rev=REV)
	output:
		#directory(expand("{out}/fastq/{target}/1/", out=OUT, target=TARGET)),
		#directory(expand("{out}/fastq/{target}/2/", out=OUT, target=TARGET)),
		out=directory("{out}/fastq/{target}/all_samples/")
	params:
		recreate_ref_folder=RECREATEREFFOLDER,
		folder=OUT,
		output="{out}/fastq/{target}/all_samples",
		forward_samples="{out}/fastq/{target}/1",
		reverse_samples="{out}/fastq/{target}/2",
		haplotype_output="{out}/haplotype_output",
		out="{out}/fastq/{target}/all_samples"
	run:
		print(input.refs)
		if params.recreate_ref_folder == True:
			if os.path.exists("ref"):
				shell("rm -r ref")
		shell("perl scripts/step1_splitSyncReadsMultiRef.pl 1 {input.refs} {params.folder} {input.pair1} {input.pair2} {input.forward} {input.rev}")
		shell("mkdir {params.output}")
		shell("mv {params.forward_samples}/*.fastq.gz {params.out} && mv {params.reverse_samples}/*.fastq.gz {params.out}")
		shell("rm -r {params.forward_samples} && rm -r {params.reverse_samples}")

rule call_haplotypes:
	input:
		"{out}/fastq/{target}/all_samples/",
	output:
		trim_filter_table="{out}/haplotype_output/{target}_trimAndFilterTable",
		results="{out}/haplotype_output/{target}_haplotypes.rds",
		reads_table="{out}/haplotype_output/{target}_trackReadsThroughPipeline.csv",
	params:
		all_samples="{out}/fastq/{target}/all_samples",
		rscript="scripts/step2_dada2.R",
		cutoff=CUTOFF,
		seed=SEED,
	script:
		"{params.rscript}"

rule censor_haplotypes:
	input:
		"{out}/haplotype_output/{target}_haplotypes.rds",
	output:
		precensored_haplotype_table="{out}/haplotype_output/{target}_haplotype_table_precensored.csv",
		snps_between_haps="{out}/haplotype_output/{target}_snps_between_haps_within_samples.fasta",
		unique_seqs="{out}/haplotype_output/{target}_uniqueSeqs.fasta",
		aligned_seqs="{out}/haplotype_output/{target}_aligned_seqs.fasta",
		final_censored="{out}/haplotype_output/{target}_uniqueSeqs_final_censored.fasta",
		final_haplotype_table="{out}/haplotype_output/{target}_haplotype_table_censored_final_version.csv",
	params:
		rscript="scripts/step3_haplotype_censoring.R",
		haplotypes="{out}/haplotype_output/{target}_haplotypes.rds",
		depth=READ_DEPTH,
		proportion=PROPORTION,
		length=HAPLOTYPE_LENGTH,
		ratio=READ_DEPTH_RATIO,
	script:
		"{params.rscript}"