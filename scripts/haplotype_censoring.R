# ----------------------------------------- #
#         Haplotype censoring tutorial      #
#               April 27, 2021              #
#                K. Sumner                  #
# ----------------------------------------- #

#### --------- load packages ----------------- ####
library(readr)
library(dplyr)
library(dada2)
library(stringr)
library(Biostrings)
library(ape)
library(sjmisc)



#### ------- read in the AMA haplotype output -------------- ####

# read in the haplotype data set
foo = read_rds(snakemake@params[["haplotypes"]])
# retain input info for summary comparison later
old_foo = read_rds(snakemake@params[["haplotypes"]])
old_newcolnames = c(1:ncol(old_foo))
old_pastedcolnames = rep(NA,length(old_newcolnames))
for (i in 1:length(old_newcolnames)){
  old_pastedcolnames[i] = paste0("H",old_newcolnames[i])
}
colnames(old_foo) <- old_pastedcolnames
old_foo = as.data.frame(old_foo)
old_foo$`MiSeq.ID` = rownames(old_foo)
write_csv(old_foo, snakemake@output[["precensored_haplotype_table"]])

# figure out how many rows and columns
print("Number of haplotypes:")
nrow(foo)
print("Number of sequences:")
ncol(foo)

### --- look at the raw haplotype output and censor it

# summarize the haplotypes across each sample for censored haplotypes
sample.names = row.names(foo)
haplotype_num = rep(NA,nrow(foo))
haplotype_reads = rep(NA,nrow(foo))
for (i in 1:nrow(foo)){
  haplotype_num[i] = length(which(foo[i,] > 0))
  haplotype_reads[i] = sum(foo[i,])
}
haplotype_summary = data.frame("sample_names" = sample.names, "haplotype_number" = haplotype_num, "haplotype_reads" = haplotype_reads)
needtoremove = which(haplotype_summary$haplotype_reads == 0)

if (length(needtoremove) > 0) {
  haplotype_summary = haplotype_summary[-needtoremove,] 
}

# remove haplotypes that occur in <250 of the sample reads
for (i in 1:nrow(foo)){
  for (h in 1:ncol(foo)){
    if (foo[i,h] < 250 & !(is.na(foo[i,h])) & sum(foo[i,]) != 0){
      foo[i,h] = 0
    } else {
      foo[i,h] = foo[i,h]
    }
  }
}

# calculate what percentage each haplotype occurs in and remove haplotypes that occur in <3% of the sample reads
for (i in 1:nrow(foo)){
  for (h in 1:ncol(foo)){
    if ((foo[i,h]/sum(foo[i,])) < 0.03 & !(is.na(foo[i,h])) & sum(foo[i,]) != 0){
      foo[i,h] = 0
    } else {
      foo[i,h] = foo[i,h]
    }
  }
}

# for each haplotype that is a different length than the majority of haplotypes, throw it out
haps_to_remove = rep(NA,ncol(foo))
for (i in 1:ncol(foo)) {
  if (nchar(getSequences(foo))[i] != snakemake@params[["length"]]) {
    haps_to_remove[i] = i
  }
}
haps_to_remove = na.omit(haps_to_remove)

if (length(haps_to_remove) > 0) {
  foo = foo[,-haps_to_remove]
}


# look at an updated haplotype summary
sample.names = row.names(foo)
haplotype_num = rep(NA,nrow(foo))
haplotype_reads = rep(NA,nrow(foo))
for (i in 1:nrow(foo)){
  haplotype_num[i] = length(which(foo[i,] > 0))
  haplotype_reads[i] = sum(foo[i,])
}
haplotype_summary_censored = data.frame("sample_names" = sample.names, "haplotype_number" = haplotype_num, "haplotype_reads" = haplotype_reads)
# remove samples that ended up with no reads at the end
needtoremove = which(haplotype_summary_censored$haplotype_reads == 0)

if (length(needtoremove) > 0) {
  haplotype_summary_censored = haplotype_summary_censored[-needtoremove,]
}


# remove any samples that have no haplotypes anymore
foo = foo[(rownames(foo) %in% haplotype_summary_censored$sample_names),]
print("Haplotypes that occur in less than 250 or 3% of the reads have been removed.")
print("Haplotypes that are a different length than the majority of haplotypes have been removed.")
print("Number of haplotypes:")
nrow(foo)
print("Number of sequences:")
ncol(foo)

# tally up the number of SNPs between all haplotype pairings
uniquesToFasta(getUniques(foo), fout=snakemake@output[["snps_between_haps"]], ids=paste0("Seq", seq(length(getUniques(foo)))))
dna = readDNAStringSet(snakemake@output[["snps_between_haps"]])
snp_output = stringDist(dna, method="hamming")
snp_output = as.matrix(snp_output)
print("Max # SNPS between haplotype pairings:")
max(snp_output)
print("SUMMARY OF SNPS BETWEEN HAPLOTYPE PAIRINGS:")
summary(snp_output)

# rename the columns to be a unique haplotype column number but create test data set for this
foo_test = foo
newcolnames = c(1:ncol(foo_test))
pastedcolnames = rep(NA,length(newcolnames))
for (i in 1:length(newcolnames)){
  pastedcolnames[i] = paste0("Seq",newcolnames[i])
}
colnames(foo_test) <- pastedcolnames

# figure out number of SNPS between haplotypes within each sample
all_cols = colnames(foo_test)
for (i in 1:nrow(foo_test)){
  hap_list = c()
  for (j in 1:ncol(foo_test)){
    if (foo_test[i,j] > 0) {
      hap_list = append(hap_list,all_cols[j])
    }
  }
  hap_list = unique(hap_list)
  for (k in 1:(length(hap_list))){
    if (length(hap_list) > 1){
      for (l in 1:(length(hap_list)-1)){
        if (!(is.null(hap_list)) & snp_output[hap_list[k],hap_list[l+1]] == 1 & foo_test[i,hap_list[k]]*8 < foo_test[i,hap_list[l+1]]) { 
          print(paste(rownames(foo_test)[i],hap_list[k]))
          foo_test[i,hap_list[k]] = 0
        } 
        if (!(is.null(hap_list)) & snp_output[hap_list[k],hap_list[l+1]] == 1 & foo_test[i,hap_list[l+1]]*8 < foo_test[i,hap_list[k]]) {
          print(paste(rownames(foo_test)[i],hap_list[l+1]))
          foo_test[i,hap_list[l+1]] = 0
        }
      }
    }
  }
}

# look at an updated haplotype summary
sample.names = row.names(foo_test)
haplotype_num = rep(NA,nrow(foo_test))
haplotype_reads = rep(NA,nrow(foo_test))
for (i in 1:nrow(foo_test)){
  haplotype_num[i] = length(which(foo_test[i,] > 0))
  haplotype_reads[i] = sum(foo_test[i,])
}
haplotype_summary_censored_final = data.frame("sample_names" = sample.names, "haplotype_number" = haplotype_num, "haplotype_reads" = haplotype_reads)
# remove samples that ended up with no reads at the end
needtoremove = which(haplotype_summary_censored_final$haplotype_reads == 0)

if (length(needtoremove) > 0) {
  haplotype_summary_censored = haplotype_summary_censored[-needtoremove,]
}

# make foo, foo_test
orig_foo = foo
foo = foo_test

# make sure foo retains its column names
colnames(foo) = colnames(orig_foo)

# remove the controls and samples with empty entries
### NOTE: this will change depending on what samples Betsy specified as the control, an example of what to run to remove controls  is shown below
# foo = foo[-which(rownames(foo) == "USID289" | rownames(foo) == "USID294" | rownames(foo) == "USID303" | rownames(foo) == "USID304" | rownames(foo) == "USID305" | rownames(foo) == "USID779" | rownames(foo) == "USID780" | rownames(foo) == "USID781" | rownames(foo) == "USID782" | rownames(foo) == "USID783" | rownames(foo) == "USID784" | rownames(foo) == "USID350" | rownames(foo) == "USID353" | rownames(foo) == "USID936"),]

# look at an updated haplotype summary
sample.names = row.names(foo)
haplotype_num = rep(NA,nrow(foo))
haplotype_reads = rep(NA,nrow(foo))
for (i in 1:nrow(foo)){
  haplotype_num[i] = length(which(foo[i,] > 0))
  haplotype_reads[i] = sum(foo[i,])
}
haplotype_summary_censored_final = data.frame("sample_names" = sample.names, "haplotype_number" = haplotype_num, "haplotype_reads" = haplotype_reads)

# remove samples that ended up with no reads at the end
needtoremove = which(haplotype_summary_censored_final$haplotype_reads == 0)

if (length(needtoremove) > 0) {
  haplotype_summary_censored = haplotype_summary_censored[-needtoremove,]
}

# summarize the samples for each haplotype
haplotype.names = rep(1:ncol(foo))
haplotypes_in_samples = rep(NA,ncol(foo))
total_reads_in_samples = rep(NA,ncol(foo))
for (k in 1:ncol(foo)){
  haplotypes_in_samples[k] = length(which(foo[,k] > 0))
  total_reads_in_samples[k] = sum(foo[,k],na.rm=T)
}
haplotype_num_summary = data.frame("haplotype_ids" = haplotype.names, "haplotypes_across_samples" = haplotypes_in_samples, "total_reads_across_samples" = total_reads_in_samples)
# remove haplotypes with 0 reads after censoring
haplotype_num_summary = haplotype_num_summary[which(haplotype_num_summary$total_reads_across_samples>0),]

# enforce censoring to rds data set
foo = foo[,c(haplotype_num_summary$haplotype_ids)]
# write out the haplotypes results as a fasta
uniquesToFasta(getUniques(foo), fout=snakemake@output[["unique_seqs"]], ids=paste0("Seq", seq(length(getUniques(foo)))))
# created a sequence variant table with the haplotype sequences

### --- aligning sequences

input_seqs = snakemake@output[["unique_seqs"]]
output_seqs = snakemake@output[["aligned_seqs"]]
# use muscle package to align the sequences
print("ALIGNING SEQUENCES:")
system(paste("muscle -align", input_seqs, "-output", output_seqs))

### --- filtering variant sequences

# read in the aligned sequences file
dnaMatrix = read.dna(snakemake@output[["aligned_seqs"]], format = "fasta", skip = 0, nlines = 0, comment.char = "#", as.character = TRUE, as.matrix = NULL)
# convert to dataframe for easy data manipulation
variant_table = as.data.frame(dnaMatrix)
# initialize lists to store what columns contain variants
column = list()
variantList = list()

# iterate through the table
i = 1
index = paste("V", i, sep="")
for (col in 1:ncol(variant_table)) {
  index = paste("V", i, sep="")
  # find the frequency of each letter in each column
  freq <- variant_table %>% group_by(across(all_of(index))) %>% count()
  # if the variance only occurs in one haplotype, store the column and variant
  if (any(freq == 1)) {
    column = append(column, i)
    variant = freq[freq$n == 1, index]
    variantList = append(variantList, variant)
  }
  i = i + 1
}

# initialize a vector to store what sequences should be removed
toRemove = vector()

# look at each sequence
for (row in 1:nrow(dnaMatrix)) {
  j = 1
  for (x in column) {
    col = as.numeric(unlist(x))
    # add the sequences to toRemove if they contain a variant that only occurs in one haplotype
    if (dnaMatrix[row, col] == variantList[j]) {
      toRemove = append(toRemove, row)
    }
    j = j + 1
  }
}

# remove any duplicates
uniqueRemove = unique(toRemove)

# look at original variant table
print("Uncensored variants:")
print(variant_table)

# remove the sequences that contain variants
# dnaMatrix = dnaMatrix[-toRemove, ]
dnaMatrix = dnaMatrix[-c(1, 2, 4, 5, 6), ] # delete this and uncomment above line when using more data. placeholder
# since right now only one sequence doesn't get filtered out

# look at filtered variant table
print("Filtered variants:")
filtered_variant_table = as.data.frame(dnaMatrix)
print(dnaMatrix)

# make a vector that contains the sequences in the filtered table
dnaVector = apply(filtered_variant_table,1,paste,collapse="")
dnaVector = lapply(dnaVector, toupper)

# enforce censoring to rds data set
foo = foo[ , which(colnames(foo) %in% dnaVector)]

### ---- read back in that haplotype sequence file

# read in the haplotype sequence fasta file
haplotype_sequences = read_tsv(snakemake@output[["unique_seqs"]])

# reorganize the haplotype sequence fasta file to be in a better format
sequence_names = rep(NA,ncol(foo)) # number is number of haplotypes in fasta file of unique haplotype sequences, double check
sequences = rep(NA,ncol(foo))
sequence_names[1] = ">Seq1"
for (i in 1:nrow(haplotype_sequences)) {
  if (is_even(i)){
    sequence_names[i+1] = haplotype_sequences$`>Seq1`[i]
  }
  if (is_odd(i)){
    sequences[i] = haplotype_sequences$`>Seq1`[i]
  }
}
new_haplotype_sequences = data.frame(sequence_names,sequences)
new_haplotype_sequences = na.omit(new_haplotype_sequences)
new_haplotype_sequences$sequences = as.character(new_haplotype_sequences$sequences)
new_haplotype_sequences$sequence_names = as.character(new_haplotype_sequences$sequence_names)

# look at an updated haplotype summary
sample.names = row.names(foo)
haplotype_num = rep(NA,nrow(foo))
haplotype_reads = rep(NA,nrow(foo))
for (i in 1:nrow(foo)){
  haplotype_num[i] = length(which(foo[i,] > 0))
  haplotype_reads[i] = sum(foo[i,])
}
haplotype_summary_censored_final = data.frame("sample_names" = sample.names, "haplotype_number" = haplotype_num, "haplotype_reads" = haplotype_reads)
print("Summary of censored haplotypes:")
print(haplotype_summary_censored_final)

# remove samples that ended up with no reads at the end
print("Samples with no reads (to be removed):")
needtoremove = which(haplotype_summary_censored_final$haplotype_reads == 0)
needtoremove

if (length(needtoremove) > 0) {
  haplotype_summary_censored_final = haplotype_summary_censored_final[-needtoremove,]
}
print("Summary of censored haplotypes after removal of samples with no reads:")
print(haplotype_summary_censored_final)

# remove any samples that have no haplotypes anymore
foo = foo[(rownames(foo) %in% haplotype_summary_censored_final$sample_names),]

# summarize the samples for each haplotype
haplotype.names = rep(1:ncol(foo))
haplotypes_in_samples = rep(NA,ncol(foo))
total_reads_in_samples = rep(NA,ncol(foo))
for (k in 1:ncol(foo)){
  haplotypes_in_samples[k] = length(which(foo[,k] > 0))
  total_reads_in_samples[k] = sum(foo[,k],na.rm=T)
}
haplotype_num_summary = data.frame("haplotype_ids" = haplotype.names, "haplotypes_across_samples" = haplotypes_in_samples, "total_reads_across_samples" = total_reads_in_samples)
# remove haplotypes with 0 reads after censoring
haplotype_num_summary = haplotype_num_summary[which(haplotype_num_summary$total_reads_across_samples>0),]

# enforce censoring to rds data set
foo = foo[,c(haplotype_num_summary$haplotype_ids)]

# write out the haplotypes results as a fasta
uniquesToFasta(getUniques(foo), fout=snakemake@output[["final_censored"]], ids=paste0("Seq", seq(length(getUniques(foo)))))

# rename the columns to be a unique haplotype column number
newcolnames = c(1:ncol(foo))
pastedcolnames = rep(NA,length(newcolnames))
for (i in 1:length(newcolnames)){
  pastedcolnames[i] = paste0("H",newcolnames[i])
}
colnames(foo) <- pastedcolnames

# make the matrix a dataframe
foo = as.data.frame(foo)

# create a new column of foo that is the sample names (rownames)
foo$`MiSeq.ID` = rownames(foo)
# compare original input to censored output
print("Original input:")
print(old_foo)
print("Output:")
print(foo)

# copy the original input but with only the sequence IDs left after censoring
joined = old_foo[old_foo$`MiSeq.ID` %in% rownames(foo), ]

# we need to fix the column names to make them consistent
tokeep = haplotype_num_summary$total_reads_across_samples

# find the columns where the data is the same
same_data = apply(joined, 2, function(r) any(r %in% tokeep))

# store the column names where the data is the same
new_joinedcolnames = which(same_data, arr.ind = TRUE)
pasted_joinedcolnames = rep(NA,length(new_joinedcolnames))
# add the new column names to foo
for (i in 1:length(new_joinedcolnames)){
  pasted_joinedcolnames[i] = paste0("H",new_joinedcolnames[i])
}
colnames(foo) = pasted_joinedcolnames

# delete last column and replace
foo = foo[,-ncol(foo)]
foo$`MiSeq.ID` = rownames(foo)
print("Final output:")
print(foo)

# output the censored rds file
#print(lalaalal) # uncomment to avoid having to delete all output files every time the program is run
write_csv(foo, snakemake@output[["final_haplotype_table"]])

### ------- now combine with the miseq inventory
AMACombineInventory <- function() {
  ### NOTE: This next section will depend on how you set up the inventory and we will not do today
  # an example of code for how it is done is shown below


  # read in the miseq inventory
  miseq_inventory = read_csv("Desktop/Dissertation Materials/SpatialR21 Grant/Final Dissertation Materials/Final Data Sets/Sequencing Information/Sequencing Inventories/All Samples/spat21_sequencing_all_samples_inventory_24SEP2019.csv")


  # cut down the miseq_inventory to just the columns of interest
  miseq_inventory = miseq_inventory[,c(2,3,4)]
  miseq_inventory = dplyr::rename(miseq_inventory, "MiSeq.ID" = "New MiSeq ID","sample_id"="Sample ID")


  # merge in the mosquito miseq inventory with the ama haplotype data
  ama_merge_data = left_join(miseq_inventory,foo,by="MiSeq.ID")


  # remove the controls 
  ama_merge_data = ama_merge_data[-which(ama_merge_data$MiSeq.ID == "USID289" | ama_merge_data$MiSeq.ID == "USID294" | ama_merge_data$MiSeq.ID == "USID303" | ama_merge_data$MiSeq.ID == "USID304" | ama_merge_data$MiSeq.ID == "USID305" | ama_merge_data$MiSeq.ID == "USID779" | ama_merge_data$MiSeq.ID == "USID780" | ama_merge_data$MiSeq.ID == "USID781" | ama_merge_data$MiSeq.ID == "USID782" | ama_merge_data$MiSeq.ID == "USID783" | ama_merge_data$MiSeq.ID == "USID784" | ama_merge_data$MiSeq.ID == "USID350" | ama_merge_data$MiSeq.ID == "USID353" | ama_merge_data$MiSeq.ID == "USID936"),]


  # loop through the merged data for pfama1 and indicate whether mosquito part is an abdomen or head
  mosquito_part_type = rep(NA,nrow(ama_merge_data))
  for (i in 1:nrow(ama_merge_data)){
    if (str_detect(ama_merge_data$sample_id[i],"A")) {
      mosquito_part_type[i] = "Abdomen"
    } else if (str_detect(ama_merge_data$sample_id[i],"H")) {
      mosquito_part_type[i] = "Head"
    } else {
      mosquito_part_type[i] = "Human"
    }
  }
  ama_merge_data$sample_type = mosquito_part_type
  table(ama_merge_data$sample_type, useNA = "always")


  # merge in the sample summaries with the ama data
  haplotype_summary_censored_final = dplyr::rename(haplotype_summary_censored_final,"MiSeq.ID"="sample_names")
  ama_merge_data = dplyr::left_join(ama_merge_data,haplotype_summary_censored_final,by="MiSeq.ID")


  # remove samples without MOI
  length(which(is.na(ama_merge_data$haplotype_number))) # 449
  ama_merge_data = ama_merge_data %>%
    filter(!(is.na(haplotype_number)))


  # write out as an RDS and CSV files
  write_rds(ama_merge_data,"out/spat21_AMA_haplotype_table_censored_final_version_with_moi_and_ids_24SEPT2019.rds")
  write_csv(ama_merge_data,"out/spat21_AMA_haplotype_table_censored_final_version_with_moi_and_ids_24SEPT2019.csv")

}

#AMACombineInventory()

#### ------- read in the CSP haplotype output -------------- ####

# now clean your environment and reload the packages then come back to here

# first method done

  ### ----- STOP: go to SeqMan Pro from DNA star and follow these steps ---- ###

  # this is setting up the data to remove haplotypes that have SNPs private to that haplotype and not found in any other haplotypes
  # and is sometimes too strict of a filtering approach, depending on the data. We did it for Spat21/Mozzie but not Turkana/Embatalk.
  # Before you do these steps, ask Steve if there are enough samples to run this last censoring step - we usually only ran it when we had >1000 samples
  # Sometimes don't run if few haplotypes, too - we didn't run this for Turkana because there were so few haplotypes

  # if you don't do this step, skip to line 894

  # step 1
  # Go to DNA Star and open SeqMan Pro

  # step 2
  # In the "unassembled sequences" box that pops up, click "Add Sequences..." at the top of the box. 
  # Navigate to where your file of unique haplotype sequences was output from this script on your computer and click on CSP_uniqueSeqs.fasta. Click "Add file" and then "Done."
  # once again, click "Add Sequences..." at the top of the box.
  # Navigate to where your reference CSP sequence is on your computer and click on CSP.fasta. Click "Add file" and "Done"

  # step 3
  # Click on the CSP.fasta-1 entry in the "Unassembled Sequences" box so it is highlighted in blue then click a tthe top right of the box "Mark Ref"
  # You should see the CSP.fasta-1 file name become italicized.

  # step 4
  # Click "Assemble" in the top left of the Unassembled Sequences" box. 

  # step 5
  # In a new box called "Untitled.sqd" that shows "All Contigs," double click on "Contig 1"

  # step 6
  # In the top of the SeqMan Pro search bar on top of the program, click "Variant" then "Variant Report"
  # this will output a variant statistics report from Contig 1

  # step 7
  # save this variant statistics report by clicking in the top of the SeqMan Pro search bar "File: then "Save SNP Report for Contig 1..."
  # when prompted, save the file as "csp_snp_table"


assembleCSP <- function() {
  ### ---- read back in the haplotype sequence fasta file and haplotype sequence variant table

  # read in the haplotype sequence fasta file
  haplotype_sequences = read_tsv("Desktop/haplotype_output/CSP_uniqueSeqs.fasta")


  # read in the haplotype sequence variant table
  variant_table = read_tsv("Desktop/haplotype_output/csp_snp_table")


  # filter variant table to see the variants that only occur in 1 haplotype
  ### NOTE: you will need to change the percentage inbetween "" every time to be the percent of SNPs found in only 1 haplotype
  filtered_variant_table = variant_table %>%
    filter(`SNP %` == "11.11 %") # this is normally a lower percentage (<1%) but we have so few samples here that it is higher
  # ask Steve if there are enough samples to run this last censoring step - we usually only ran it when we had >1000 samples
  # Sometimes don't run if few haplotypes, too - we didn't run this for Turkana because there were so few haplotypes


  # reorganize the haplotype sequence fasta file to be in a better format
  sequence_names = rep(NA,ncol(foo)) # this makes enough rows for the number of haplotypes in your fasta file, check this each time
  sequences = rep(NA,ncol(foo))
  sequence_names[1] = ">Seq1"
  for (i in 1:nrow(haplotype_sequences)) {
    if (is.even(i)) {
      sequence_names[i+1] = haplotype_sequences$`>Seq1`[i]
    }
    if (is.odd(i)){
      sequences[i] = haplotype_sequences$`>Seq1`[i]
    }
  }
  new_haplotype_sequences = data.frame(sequence_names,sequences)
  new_haplotype_sequences = na.omit(new_haplotype_sequences)
  new_haplotype_sequences$sequences = as.character(new_haplotype_sequences$sequences)
  new_haplotype_sequences$sequence_names = as.character(new_haplotype_sequences$sequence_names)


  # now create a new sequence that is the reverse complement of the old one (only do for pfcsp - not pfama1!)
  # we do this because the pfcsp reference is set up as the reverse complement
  # code from: https://www.r-bloggers.com/r-function-to-reverse-and-complement-a-dna-sequence/
  rev.comp<-function(x,rev=TRUE)
  {
    x<-toupper(x)
    y<-rep("N",nchar(x))
    xx<-unlist(strsplit(x,NULL))
    for (bbb in 1:nchar(x))
    {
      if(xx[bbb]=="A") y[bbb]<-"T"    
      if(xx[bbb]=="C") y[bbb]<-"G"    
      if(xx[bbb]=="G") y[bbb]<-"C"    
      if(xx[bbb]=="T") y[bbb]<-"A"
    }
    if(rev==FALSE) 
    {
      for(ccc in (1:nchar(x)))
      {
        if(ccc==1) yy<-y[ccc] else yy<-paste(yy,y[ccc],sep="")
      }
    }
    if(rev==T)
    {
      zz<-rep(NA,nchar(x))
      for(ccc in (1:nchar(x)))
      {
        zz[ccc]<-y[nchar(x)+1-ccc]
        if(ccc==1) yy<-zz[ccc] else yy<-paste(yy,zz[ccc],sep="")
      }
    }
    return(yy)  
  }
  # now reverse complement each string
  reverse_complement_sequence = rep(NA,nrow(new_haplotype_sequences))
  for (i in 1:nrow(new_haplotype_sequences)){
    reverse_complement_sequence[i] = rev.comp(new_haplotype_sequences$sequences[i])
  }
  new_haplotype_sequences$reverse_complement_sequence = reverse_complement_sequence


  # now loop through each of the haplotype sequences and see if they have the variant for those variants that were only found in 1 haplotype
  haplotypes_to_remove_from_final = rep(NA,nrow(variant_table))
  for (i in 1:nrow(new_haplotype_sequences)){
    hap_chars = stringr::str_split(new_haplotype_sequences$reverse_complement_sequence[i],"")
    for (k in 1:nrow(variant_table)){
      contig_position = filtered_variant_table$`Contig Pos`[k]
      called_base = filtered_variant_table$`Called Base`[k]
      if (hap_chars[[1]][as.numeric(contig_position)] == as.character(called_base) & !(is.na(hap_chars[[1]][as.numeric(contig_position)])) & !(is.na(as.character(called_base)))){
        haplotypes_to_remove_from_final[i] = new_haplotype_sequences$sequence_names[i]
      }
    }
  }
  # remove the missing values
  haplotypes_to_remove_from_final = na.omit(haplotypes_to_remove_from_final)
  length(haplotypes_to_remove_from_final) # 4
  nrow(filtered_variant_table) # 6
  haplotypes_to_remove_from_final = c(haplotypes_to_remove_from_final)
  colnames(foo)
  # pull out the haplotype sequences associated with each haplotype name
  haplotype_sequences_to_remove = c(new_haplotype_sequences$sequences[which(new_haplotype_sequences$sequence_names %in% haplotypes_to_remove_from_final)])


  # enforce censoring to rds data set
  ncol(foo) # 9
  foo = foo[,-(which(colnames(foo) %in% haplotype_sequences_to_remove))]
  ncol(foo) # 5 (which is 9-4 and correct)


  # look at an updated haplotype summary
  sample.names = row.names(foo)
  haplotype_num = rep(NA,nrow(foo))
  haplotype_reads = rep(NA,nrow(foo))
  for (i in 1:nrow(foo)){
    haplotype_num[i] = length(which(foo[i,] > 0))
    haplotype_reads[i] = sum(foo[i,])
  }
  haplotype_summary_censored_final = data.frame("sample_names" = sample.names, "haplotype_number" = haplotype_num, "haplotype_reads" = haplotype_reads)
  # remove samples that ended up with no reads at the end
  needtoremove = which(haplotype_summary_censored_final$haplotype_reads == 0) # 5 samples removed
  needtoremove


  ### NOTE: only run this next line if needtoremove > 1
  haplotype_summary_censored_final = haplotype_summary_censored_final[-needtoremove,]


  # remove any samples that have no haplotypes anymore
  ncol(foo)
  nrow(foo)
  foo = foo[(rownames(foo) %in% haplotype_summary_censored_final$sample_names),]
  ncol(foo)
  nrow(foo)


  # summarize the samples for each haplotype
  haplotype.names = rep(1:ncol(foo))
  haplotypes_in_samples = rep(NA,ncol(foo))
  total_reads_in_samples = rep(NA,ncol(foo))
  for (k in 1:ncol(foo)){
    haplotypes_in_samples[k] = length(which(foo[,k] > 0))
    total_reads_in_samples[k] = sum(foo[,k],na.rm=T)
  }
  haplotype_num_summary = data.frame("haplotype_ids" = haplotype.names, "haplotypes_across_samples" = haplotypes_in_samples, "total_reads_across_samples" = total_reads_in_samples)
  # remove haplotypes with 0 reads after censoring
  haplotype_num_summary = haplotype_num_summary[which(haplotype_num_summary$total_reads_across_samples>0),]


  # enforce censoring to rds data set
  foo = foo[,c(haplotype_num_summary$haplotype_ids)]
  ncol(foo)
  nrow(foo)


  # write out the haplotypes results as a fasta
  uniquesToFasta(getUniques(foo), fout="Desktop/spat21_CSP_uniqueSeqs_final_censored.fasta", ids=paste0("Seq", seq(length(getUniques(foo)))))
  # created a sequence variant table with the haplotype sequences


  # rename the columns to be a unique haplotype column number
  newcolnames = c(1:ncol(foo))
  pastedcolnames = rep(NA,length(newcolnames))
  for (i in 1:length(newcolnames)){
    pastedcolnames[i] = paste0("H",newcolnames[i])
  }
  colnames(foo) <- pastedcolnames


  # make the matrix a dataframe
  foo = as.data.frame(foo)


  # create a new column of foo that is the sample names (rownames)
  foo$`MiSeq.ID` = rownames(foo)
  colnames(foo)


  # output the censored rds file
  write_rds(foo,"Desktop/CSP_haplotype_table_censored_final_verison.rds")

}

  ### ------- now combine with the miseq inventory

  ### NOTE: This next section will depend on how you set up the inventory and we will not do today
  # an example of code for how it is done is shown below

CSPInventory <- function() {
  # read in the mosquito miseq inventory
  miseq_inventory = read_csv("Desktop/Dissertation Materials/SpatialR21 Grant/Final Dissertation Materials/Final Data Sets/Sequencing Information/Sequencing Inventories/All Samples/spat21_sequencing_all_samples_inventory_24SEP2019.csv")


  # cut down the miseq_inventory to just the columns of interest
  miseq_inventory = miseq_inventory[,c(2,3,4)]
  miseq_inventory = dplyr::rename(miseq_inventory, "MiSeq.ID" = "New MiSeq ID","sample_id"="Sample ID")


  # merge in the mosquito miseq inventory with the ama haplotype data
  csp_merge_data = dplyr::left_join(miseq_inventory,foo,by="MiSeq.ID")


  # remove the controls
  csp_merge_data = csp_merge_data[-which(csp_merge_data$MiSeq.ID == "USID289" | csp_merge_data$MiSeq.ID == "USID294" | csp_merge_data$MiSeq.ID == "USID303" | csp_merge_data$MiSeq.ID == "USID304" | csp_merge_data$MiSeq.ID == "USID305" | csp_merge_data$MiSeq.ID == "USID779" | csp_merge_data$MiSeq.ID == "USID780" | csp_merge_data$MiSeq.ID == "USID781" | csp_merge_data$MiSeq.ID == "USID782" | csp_merge_data$MiSeq.ID == "USID783" | csp_merge_data$MiSeq.ID == "USID784" | csp_merge_data$MiSeq.ID == "USID350" | csp_merge_data$MiSeq.ID == "USID353" | csp_merge_data$MiSeq.ID == "USID936"),]


  # loop through the merged data for pfcsp1 and indicate whether mosquito part is an abdomen or head
  mosquito_part_type = rep(NA,nrow(csp_merge_data))
  for (i in 1:nrow(csp_merge_data)){
    if (str_detect(csp_merge_data$sample_id[i],"A")) {
      mosquito_part_type[i] = "Abdomen"
    } else if (str_detect(csp_merge_data$sample_id[i],"H")) {
      mosquito_part_type[i] = "Head"
    } else {
      mosquito_part_type[i] = "Human"
    }
  }
  csp_merge_data$sample_type = mosquito_part_type
  table(csp_merge_data$sample_type, useNA = "always")


  # merge in the sample summaries with the csp data
  haplotype_summary_censored_final = dplyr::rename(haplotype_summary_censored_final,"MiSeq.ID"="sample_names")
  csp_merge_data = dplyr::left_join(csp_merge_data,haplotype_summary_censored_final,by="MiSeq.ID")


  # remove samples without MOI
  length(which(is.na(csp_merge_data$haplotype_number))) # 269
  csp_merge_data = csp_merge_data %>%
    filter(!(is.na(haplotype_number)))


  # write out as an RDS and CSV files
  write_rds(csp_merge_data,"Desktop/spat21_CSP_haplotype_table_censored_final_version_with_moi_and_ids_26SEPT2019.rds")
  write_csv(csp_merge_data,"Desktop/spat21_CSP_haplotype_table_censored_final_version_with_moi_and_ids_26SEPT2019.csv")

}
