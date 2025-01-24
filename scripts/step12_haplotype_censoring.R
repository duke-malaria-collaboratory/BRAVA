###################################
#                                 #
#      Haplotype censoring        #
#    C Markwalter and Z Lapp      #
#           June 2023             #
#                                 #
###################################

#### Packages
library(msa)
library(ape)
library(reshape2)
library(tidyverse)

#### Read data in and set parameters
precensored <- read_csv(snakemake@input[[1]]) %>%
  replace(is.na(.), 0)

ratio_thresh <- 1/snakemake@params$ratio
hap_length <- read_csv(snakemake@input[[2]]) %>% select(snakemake@wildcards$target) %>% unlist()
read_depth <- snakemake@params$depth
read_prop <- snakemake@params$proportion

#### Read ratio
# for any haplotype (k) that is only 1 nucleotide different from another (j) in the sample and also is present at read ratio of k/j < read_ratio, then merge the reads for k into j.

# make dna stringset
dna_ss <- DNAStringSet(colnames(precensored %>% select(-seq_id)))

# set names of haplotypes as sequences
names(dna_ss) <- colnames(precensored %>% select(-seq_id))

# set penalties for mismatches
sigma <- nucleotideSubstitutionMatrix(match = 2, mismatch = -1, baseOnly = TRUE)

# multiple sequence alignment
seqs_msa <- msa(dna_ss, substitutionMatrix = sigma, gapOpening = -8, gapExtension = -5)

# get distances between sequences
dists <- dist.dna(msaConvert(seqs_msa, 'ape::DNAbin'), 
                  model = 'N', pairwise.deletion = TRUE, as.matrix = TRUE) %>% 
  melt() %>% 
  filter(Var1 != Var2) %>%
  mutate(len1 = nchar(as.character(Var1)),
         len2 = nchar(as.character(Var2))) %>%
  filter(len1 == len2)

# within a sample, check for haplotype pairs with distances of only 1 nucleotide
precensored_long <- precensored %>%
  pivot_longer(-c(seq_id), names_to = "hap", values_to = "reads") %>%
  filter(reads > 0) %>%
  group_by(seq_id) %>%
  mutate(length = nchar(hap))

precensored_long <- precensored_long %>%
  rowwise() %>%
  mutate(min_dist = ifelse(sum(dists$Var1 == hap) == 0, NA, dists %>%
                             filter(Var1 == hap & Var2 %in% precensored_long$hap[precensored_long$seq_id == seq_id]) %>%
                             slice_min(value) %>%
                             pull(value) %>%
                             unique()),
         min_dist_hap = ifelse(is.na(min_dist), NA,
                               dists %>% 
                                 filter(Var1 == hap & Var2 %in% precensored_long$hap[precensored_long$seq_id == seq_id]) %>%
                                 slice_min(value) %>%
                                 mutate(Var2 = as.character(Var2)) %>%
                                 pull(Var2) %>%
                                 unique()),
         read_ratio = ifelse(is.na(min_dist) | min_dist != 1, NA,
                             reads/precensored_long$reads[precensored_long$seq_id == seq_id & precensored_long$hap == min_dist_hap])) %>%
  ungroup() %>%
  mutate(read_ratio = ifelse(read_ratio > 1, NA, read_ratio))

low_ratio <- precensored_long %>%
  filter(read_ratio < ratio_thresh)

high_ratio_merged <- precensored_long %>%
  filter(seq_id %in% low_ratio$seq_id) %>%
  left_join(low_ratio %>%
              select(-hap, -length, -min_dist) %>%
              rename(hap = min_dist_hap,
                     ratio_low = read_ratio,
                     reads_low = reads)) %>%
  filter(!is.na(reads_low)) %>%
  group_by(seq_id, hap, reads, length) %>%
  summarise(reads = reads + sum(reads_low)) %>%
  distinct()

post_read_ratio <- precensored_long %>%
  filter(read_ratio >= ratio_thresh | is.na(read_ratio)) %>%
  select(-min_dist_hap, -min_dist, -read_ratio) %>%
  mutate(note = NA) %>%
  bind_rows(high_ratio_merged %>%
              mutate(note = "merged with low read ratio")) %>%
  group_by(seq_id, hap) %>%
  mutate(n = n()) %>%
  filter(n == 1 | (n == 2 & !is.na(note))) %>%
  ungroup()


#### Filter only to the correct length, read depth, and within-sample read proportion

censored <- post_read_ratio %>%
  filter(length == hap_length, reads > read_depth) %>% #filter out haps with incorrect length and/or low read depths
  group_by(seq_id) %>%
  mutate(prop_in_sample = reads/sum(reads)) %>%
  filter(prop_in_sample > read_prop) %>% #filter out low prop_in_sample
  select(-n, -prop_in_sample)


#### add notes to precensored haps for export

precensored_out <- precensored_long %>%
  mutate(wrong_length = ifelse(length != hap_length, 1, 0),
         low_reads = ifelse(reads < read_depth, 1, 0),
         read_ratio_low = ifelse(read_ratio < ratio_thresh, 1, 0),
         read_prop = NA) %>%
  bind_rows(precensored_long %>%
              filter(!paste0(seq_id, "_", hap) %in% paste0(censored$seq_id, "_", censored$hap),
                     length == hap_length,
                     reads > read_depth,
                     (read_ratio >= ratio_thresh | is.na(read_ratio))) %>%
              mutate(wrong_length = 0,
                     low_reads = 0,
                     read_ratio_low = 0,
                     read_prop = 1)
              ) %>%
  group_by(seq_id, hap) %>%
  mutate(n = n()) %>%
  ungroup() %>%
  filter(n == 1 | (n == 2 & read_prop == 1)) %>%
  mutate(reason_filtered = case_when(read_ratio_low == 1 ~ "ratio vs hap with dist of 1 < thresh",
                                     wrong_length == 1 & low_reads == 1 ~ "wrong length, low reads",
                                     wrong_length == 1 ~ "wrong length",
                                     low_reads == 1 ~ "low reads",
                                     read_prop == 1 ~ "proportion of reads in sample low")) %>%
  select(-n)


#### Write files
write_csv(precensored_out, snakemake@output[[1]])
write_csv(censored, snakemake@output[[2]])

