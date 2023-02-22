targets <- snakemake@params[["targets"]]
out <- snakemake@params[["out"]]

columns <- c("Target","Sample","POS","REF","ALT","ALT_DEPTH","DEPTH","ALT_FREQ") 
df <- data.frame(matrix(nrow = 0, ncol = length(columns))) 
colnames(df) <- columns

for (target in targets) {
    path <- paste0(target, "/", out, "/DR_mutations.csv")
    target_table = read.csv(path)
    target_table <- target_table[, -1]
    df <- rbind(df,target_table)
}
print(df)
write.csv(df, snakemake@output[[1]])