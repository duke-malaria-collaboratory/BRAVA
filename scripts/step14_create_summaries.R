library(tidyverse)
library(readxl)

# trim summary work
trim_summary_names <- read_delim(snakemake@params[["trim_summary_names"]], 
                                 delim = ' ', col_names = 'sample') %>% 
  mutate(sample = gsub('.summary', '', sample))

long_trim_summary <- read_delim(snakemake@params[["trim_summaries"]], 
                           col_names = c('read_type', 'read_ct')) %>% 
  mutate(read_ct = as.numeric(read_ct),
         type = ifelse(grepl('Percent', read_type), 'Percent', 'Count'),
         read_type = gsub(' Percent', 's', read_type)) %>% 
  bind_cols(trim_summary_names %>% slice(rep(1:n(), each = 9))) %>% 
  filter(type == 'Count') %>%
  select(-type) 

print("Long trim summary:")
print(long_trim_summary)

# read count summary work
pre_filt_read_cts <- read_delim(snakemake@params[["pre_filt_fastq_counts"]],
                                col_names = c('path', 'read_ct')) %>%
  mutate(sample = gsub('.*/|_..fastq.gz', '', path),
         target = gsub('.*_', '', sample),
         sample = gsub('_.*', '', sample),
         read_type = 'prefilt', 
         .before = 1) %>% 
  select(-path) %>% 
  mutate(target = tolower(target))

filt_read_cts <- read_delim(snakemake@params[["filt_fastq_counts"]],
                            col_names = c('path', 'read_ct')) %>%
  mutate(sample = gsub('.*/|_._filt.fastq.gz|_final', '', path),
         target = str_match(path, 'haplotype_output/(\\w+?)/bb')[,2],
         read_type = 'filt', 
         .before = 1) %>% 
  select(-path) %>% 
  mutate(target = tolower(target))

long_read_cts <- bind_rows(pre_filt_read_cts, filt_read_cts) %>% 
  distinct() %>%
  group_by(sample, target, read_type) %>% 
  slice_max(read_ct)

print(long_read_cts)

print("Long summary:")
long_summary <- bind_rows(long_trim_summary, long_read_cts)
print(long_summary, n=Inf)
print("Wide summary:")
wide_summary <- long_summary %>% pivot_wider(names_from = c(target, read_type), values_from = read_ct) %>%
  dplyr::rename_with(~str_remove(., "NA_"))
print(wide_summary, n=Inf)

write_csv(long_summary, snakemake@output[["long_summary"]])
write_csv(wide_summary, snakemake@output[["wide_summary"]])
