library(tidyverse)

targets <- snakemake@params[["targets"]]

columns <- c("Target","Sample","POS","REF","ALT","ALT_DEPTH","DEPTH","ALT_FREQ") 
df <- data.frame(matrix(nrow = 0, ncol = length(columns))) 
colnames(df) <- columns

for (target in targets) {
    path <- paste0(target, "/out/DR_mutations.csv")
    target_table = read.csv(path)
    target_table <- target_table[, -1]
    df <- rbind(df,target_table)
}

# dataframe of positions of interest 

pos_of_interest <- read_csv(snakemake@params[["table"]]) %>%
  mutate(Target = toupper(Target)) %>% # this is to make the targets match between this file and the vcf - MIGHT HAVE TO CHANGE THIS IF IT BREAKS FOR SOME GENES 
  filter(Target %in% targets) # filter to targets that were run through the pipeline

# vector of positions of interest
pos <- pos_of_interest %>% mutate(pos_base = paste0(Target, '_', Pos)) %>% select(pos_base)

# vcf 
vcf <- df %>%
  # only keep positions that we're interested in
  filter(paste0(Target, '_', POS) %in% pull(pos))

# bases of interest
pos_base <- bind_rows(pos_of_interest %>% 
                        mutate(pos_base = paste0(Target, '_', Pos, Ref)) %>% 
                        select(pos_base),
                      pos_of_interest %>% 
                        mutate(pos_base = paste0(Target, '_', Pos, Alt)) %>% 
                        select(pos_base))

all_dat <- bind_rows(vcf %>% 
            filter(ALT != '.') %>% 
  select(Target, Sample, POS, ALT, DEPTH, ALT_DEPTH, ALT_FREQ) %>% 
  rename(BASE = ALT, TOTAL_DEPTH = DEPTH, DEPTH = ALT_DEPTH, FREQ = ALT_FREQ),
  vcf %>% 
    filter(ALT != '.') %>% 
    mutate(REF_DEPTH = DEPTH - ALT_DEPTH,
           REF_FREQ = 1 - ALT_FREQ) %>% 
    select(Target, Sample, POS, REF, DEPTH, REF_DEPTH, REF_FREQ) %>% 
    rename(BASE = REF, TOTAL_DEPTH = DEPTH, DEPTH = REF_DEPTH, FREQ = REF_FREQ),
  vcf %>% 
    filter(ALT == '.') %>% 
    mutate(REF_DEPTH = DEPTH - ALT_DEPTH,
           REF_FREQ = 1 - ALT_FREQ) %>% 
    select(Target, Sample, POS, REF, DEPTH, REF_DEPTH, REF_FREQ) %>% 
    rename(BASE = REF, TOTAL_DEPTH = DEPTH, DEPTH = REF_DEPTH, FREQ = REF_FREQ)) #%>% 
   # mutate(pos_base = paste0(Target, '_', POS, BASE)) %>%
  # add in any positions that are not found in any samples
  #full_join(pos_base) #%>% 
  # select(-c(Target, POS, BASE)) %>% 
  # pivot_wider(names_from = pos_base, values_from = c(DEPTH, TOTAL_DEPTH, FREQ),
  #             values_fill = 0, 
  #             names_vary = 'slowest', names_glue = "{pos_base}_{.value}", names_sort = TRUE) %>% 
  #filter(!is.na(Sample))

all_dat %>% 
  write_csv(snakemake@output[["dr_depths"]])

# all_dat %>% 
#   select_if(grepl('FREQ', colnames(.))) %>% 
#   write_csv(snakemake@output[["dr_freqs"]])

unlink(snakemake@params[["DR_mutations"]])