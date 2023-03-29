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
path <- snakemake@params[["mapped_reads"]]
print("Fastq files:")
list.files(path)

# filter and trim the fastq files
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

  # performing filtering and trimming
  # first define the filenames for the filtered fastq.gz files
  filt_path <- file.path(path, "filtered") # Place filtered files in filtered/ subdirectory
  filtFs <- file.path(filt_path, paste0(sample.names, "_", snakemake@params[["q_values"]], "_F_filt.fastq.gz"))
  filtRs <- file.path(filt_path, paste0(sample.names, "_", snakemake@params[["q_values"]], "_R_filt.fastq.gz"))

  # filter the forward and reverse reads
  # remember that dada2 requires no Ns, so maxN makes sure sequences with more than maxN Ns will be discarded
  out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, maxN=0, maxEE=c(1,1), truncQ=strtoi(snakemake@params[["q_values"]]), rm.phix=FALSE,
                      compress=TRUE, multithread=FALSE, matchIDs = TRUE)
  # sum read counts and output it into a file
  sum <- colSums(out)
  read_count <- sum[['reads.out']]
  print("Number of reads:")
  print(read_count)
  q_values <- snakemake@params[["q_values"]]
  to_write <- data.frame(q_values, read_count)
  print(to_write)
  if (!dir.exists(snakemake@params[["trim_filter_out"]])) {
    dir.create(snakemake@params[["trim_filter_out"]])
    write.table(to_write, file=snakemake@params[["read_count"]], append=FALSE, col.names = FALSE, row.names = FALSE)
  } else {
    if (!file.exists(snakemake@params[["read_count"]])) {
      file.create(snakemake@params[["read_count"]])
      write.table(to_write, file=snakemake@params[["read_count"]], append=FALSE, col.names = FALSE, row.names = FALSE)
    } else {
      write.table(to_write, file=snakemake@params[["read_count"]], append=TRUE, col.names = FALSE, row.names = FALSE)
    }
  }
  # output summary of read trimming and filtering
  write.csv(out,snakemake@output[["trim_filter_table"]])
}