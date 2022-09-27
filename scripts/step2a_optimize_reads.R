library(dplyr)
require(data.table)

# read in read_counts
counts <- read.table(snakemake@params[["read_count"]])
counts <- counts %>% setNames(c("q_value", "read_count"))
# find the q value with the highest read count
max_read_count <- counts[which.max(counts$read_count),]
final_q_value <- max_read_count$q_value
final_read_count <- max_read_count$read_count
# write the highest read count to a file
write(final_read_count, file=snakemake@output[["max_read_count"]])
write(final_q_value, file=snakemake@output[["final_q_value"]])
# copy the trim and filter table for the q value with the highest read count - we will later delete all trim and filter tables except this
trim_filter_table <- file.path(snakemake@params[["target"]], snakemake@params[["out"]], "temp_haplotype_output", paste0(snakemake@params[["target"]], "_", final_q_value, "_trimAndFilterTable"))
file.copy(from=trim_filter_table, to=snakemake@output[["final_trim_filter_table"]])

# copy all the filtered fastq files for the q value with the highest read count - we will later delete all filtered fastq files except this
# create new folder for final filtered fastq files
dir.create(snakemake@params[["final_filtered"]])
path <- snakemake@params[["all_samples"]]
fnFs <- sort(list.files(path, pattern="_1.fastq.gz"))
fnRs <- sort(list.files(path, pattern="_2.fastq.gz"))
sample.names <- sapply(strsplit(fnFs, "_"), `[`, 1)
fnFs <- file.path(path, fnFs)
fnRs <- file.path(path, fnRs)
filt_path <- file.path(path, "filtered")
final_filt_path <- file.path(snakemake@params[["final_filtered"]])
filtFs <- file.path(filt_path, paste0(sample.names, "_", final_q_value, "_F_filt.fastq.gz"))
filtRs <- file.path(filt_path, paste0(sample.names, "_", final_q_value, "_R_filt.fastq.gz"))
final_filtFs <- file.path(final_filt_path, paste0(sample.names, "_final_F_filt.fastq.gz"))
final_filtRs <- file.path(final_filt_path, paste0(sample.names, "_final_R_filt.fastq.gz"))
file.copy(from=filtFs, to=final_filtFs)
file.copy(from=filtRs, to=final_filtRs)

# now i have to delete all the files