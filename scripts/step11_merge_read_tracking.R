library(dplyr)

final_q_value <- read.table(snakemake@params[["final_q_value"]])
track_reads_through_pipeline <- read.table(snakemake@params[["track_reads_through_pipeline"]], header = TRUE, sep = ",")
reads_table <- read.table(file.path(paste0(snakemake@params[["trim_filter_path"]], "_", final_q_value, "_trim_and_filter_table")), header = TRUE, sep = ",")
reads_table[, 1] <- sub("\\_.*", "", reads_table[, 1])
colnames(reads_table)[1] <- "sample"
colnames(track_reads_through_pipeline)[1] <- "sample"
print(track_reads_through_pipeline)
print(reads_table)
track_reads_through_dada2 <- full_join(reads_table, track_reads_through_pipeline, by="sample")
print("dada2:")
print(track_reads_through_dada2)
write.csv(track_reads_through_dada2, snakemake@output[["track_reads_through_dada2"]])