# ----------------------------------------- #
#       Haplotype Calling - Embatalk        #
#                Phase 1 and 2              #
#               March 31, 2020              #
#                K. Sumner                  #
# ----------------------------------------- #


# directions for this tutorial can be found at: http://benjjneb.github.io/dada2/tutorial.html


#### --------- load packages ----------------- ####

# load the dada2 library
library("dada2")

# read in the path to your folder of fastq files
path <- snakemake@params[["all_samples"]]
trim_filter_table <- read.csv(snakemake@params[["trim_filter_table"]])
print("Fastq files:")
list.files(path)

# Sort ensures forward/reverse reads are in same order
fnFs <- sort(list.files(path, pattern="_1.fastq.gz"))
fnRs <- sort(list.files(path, pattern="_2.fastq.gz"))

# check the files first
if (!any(duplicated(c(fnFs, fnRs)))) {

  # Extract sample names, assuming filenames have format: SAMPLENAME_X.fastq
  # note: the string manipulations may have to be modified, especially the extraction of sample names from the file names
  sample.names <- sapply(strsplit(fnFs, "_"), `[`, 1)
  # Specify the full path to the fnFs and fnRs
  fnFs <- file.path(path, fnFs)
  fnRs <- file.path(path, fnRs)

  # first define the filenames for the filtered fastq.gz files
  filt_path <- file.path(path, "final_filtered") # Place filtered files in filtered/ subdirectory
  filtFs <- file.path(filt_path, paste0(sample.names, "_final_F_filt.fastq.gz"))
  filtRs <- file.path(filt_path, paste0(sample.names, "_final_R_filt.fastq.gz"))

  # remove samples that had less than 50 reads after sampling
  keep <- trim_filter_table[,"reads.out"] > snakemake@params[["cutoff"]] # Or other cutoff
  filtFs <- filtFs[keep]
  filtRs <- filtRs[keep]

  # pull out sample names for the filtered reads
  sample.names <- sapply(strsplit(basename(filtFs), "_"), `[`, 1)
  sample.namesR <- sapply(strsplit(basename(filtRs), "_"), `[`, 1)
  if(!identical(sample.names, sample.namesR)) stop("Forward and reverse files do not match.")
  names(filtFs) <- sample.names
  names(filtRs) <- sample.names
  # set the random seed
  set.seed(snakemake@params[["seed"]])
  # learn the error rates for your amplicon data set
  errF <- learnErrors(filtFs, multithread=TRUE)
  errR <- learnErrors(filtRs, multithread=TRUE)

  # Sample inference and merger of paired-end reads
  sample.names <- names(filtFs)
  mergers <- vector("list", length(sample.names))
  names(mergers) <- sample.names
  for(sam in sample.names) {
    cat("Processing:", sam, "\n")
    derepF <- derepFastq(filtFs[[sam]])
    ddF <- dada(derepF, err=errF, multithread=TRUE)
    derepR <- derepFastq(filtRs[[sam]])
    ddR <- dada(derepR, err=errR, multithread=TRUE)
    merger <- mergePairs(ddF, derepF, ddR, derepR)
    mergers[[sam]] <- merger
  }
  rm(derepF); rm(derepR)

  # Construct sequence table
  seqtab <- makeSequenceTable(mergers)

  # remove the chimeras
  seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
  print("Dimensions of results after chimera removal:")
  print(dim(seqtab.nochim))
  print("Percentage of reads left after chimera removal:")
  print(sum(seqtab.nochim)/sum(seqtab))
  # note: most of your reads should remain after chimera removal
  # it is not uncommon for a majority of the sequence variants to be removed (which we observed)

  # write out the results without chimeras
  saveRDS(seqtab.nochim, snakemake@output[["results"]])

  # track reads through the pipeline
  # look at the number of reads that made it through each step in the pipeline
  ## NOTE: not working for when you have filtered reads header so took that out (have number of reads filtered in previous table)
  getN <- function(x) sum(getUniques(x))
  track <- cbind(sapply(mergers, getN), rowSums(seqtab), rowSums(seqtab.nochim))
  colnames(track) <- c("merged", "tabled", "nonchim")
  rownames(track) <- sample.names
  write.csv(track, snakemake@output[["reads_table"]])

}