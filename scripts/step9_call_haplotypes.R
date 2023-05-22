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
library(tidyverse)

# read in the path to your folder of fastq files
path <- snakemake@params[["mapped_reads"]]
trim_filter_table <- read.csv(snakemake@params[["trim_filter_table"]])
q_values <- snakemake@params[["q_values"]]
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
  filt_path <- file.path(path, "filtered") # Place filtered files in filtered/ subdirectory
  filtFs <- file.path(filt_path, paste0(sample.names, "_", q_values, "_F_filt.fastq.gz"))
  filtRs <- file.path(filt_path, paste0(sample.names, "_", q_values, "_R_filt.fastq.gz"))

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
  print(length(filtFs))
  print(length(filtRs))
# set the random seed
  set.seed(snakemake@params[["seed"]])
  # learn the error rates for your amplicon data set
  print("before errF")
  errF <- learnErrors(filtFs, multithread=TRUE)
  print("after errF")
  errR <- learnErrors(filtRs, multithread=TRUE)
  print("after errR")

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
    print("pairs:")
    print(derepF)
    print(ddF)
    print(derepR)
    print(ddR)
    merger <- mergePairs(ddF, derepF, ddR, derepR)
    print(merger)
    mergers[[sam]] <- merger
  }
  rm(derepF); rm(derepR)
  # print("mergers")
  # print(mergers)
  # Construct sequence table
  seqtab <- makeSequenceTable(mergers)
  # print("seqtab")
  # print(seqtab)
  # remove the chimeras
  seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
  print("Dimensions of results after chimera removal:")
  print(dim(seqtab.nochim))
  print("Percentage of reads left after chimera removal:")
  print(sum(seqtab.nochim)/sum(seqtab))
  # note: most of your reads should remain after chimera removal
  # it is not uncommon for a majority of the sequence variants to be removed (which we observed)

  # write out the results without chimeras
  write_csv(data.frame(seqtab.nochim)  %>% rownames_to_column(var = "seq_id"), snakemake@output[["results"]])

  # track reads through the pipeline
  # look at the number of reads that made it through each step in the pipeline
  ## NOTE: not working for when you have filtered reads header so took that out (have number of reads filtered in previous table)
  getN <- function(x) sum(getUniques(x))
  track <- cbind(sapply(mergers, getN), rowSums(seqtab), rowSums(seqtab.nochim))
  colnames(track) <- c("merged", "tabled", "nonchim")
  rownames(track) <- sample.names

  # sum read counts and output it to a file
  sum <- colSums(track)
  read_count <- sum[['nonchim']]
  print("Number of reads:")
  print(read_count)
  to_write <- data.frame(q_values, read_count)
  print(to_write)
  if (!file.exists(snakemake@params[["read_count"]])) {
    file.create(snakemake@params[["read_count"]])
    write.table(to_write, file=snakemake@params[["read_count"]], append=FALSE, col.names = FALSE, row.names = FALSE)
  } else {
      write.table(to_write, file=snakemake@params[["read_count"]], append=TRUE, col.names = FALSE, row.names = FALSE)
  }

  write.csv(track, snakemake@output[["reads_table"]])

}
