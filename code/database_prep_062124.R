
### VZV TCRs from Wen 2024

library(readr)
library(tidyverse)
library(dplyr)
library(tidyr)
library(readxl)
library(data.table)
library(stringr)
library(gdata)
library(purrr)

rm(list=ls())

# INPUT
user = 'esf'
if(user == 'esf') {
  repo_loc = '/Volumes/corey_l/esford3_kmayerbl_collab/software/wen_2024'
} else {
  stop("set repo loc and repo manually")
}

filename_main = file.path(repo_loc, 'data/wen_nc_060624.csv')

# OUTPUT
filename_out = file.path(repo_loc, 'data/wen_for_tcrdist_060724.csv')


#### DATA DICTIONARY ####
#libid - unique chain read id
#full_nt_sequence - nucleotide sequence
#v_gene - v gene
#j_gene - j gene
#junction - cdr3 amino acid sequence
#chain - alpha/beta 
#tetDonorID - HLA-matched tetramer isolation 
#donor - donor
#visit - sample collection time point 


# PREP FILE FOR CLUSTERING ANALYSIS
wen <- read_csv(filename_main)
colnames(wen)

pairs <- wen %>%
  group_by(libid) %>%
  filter(n() == 2) %>%
  filter(sum(chain == "a") < 2 & sum(chain == "b") < 2) %>%
  filter(chain == "a" | chain == "b") %>% 
  mutate(tet = str_split_i(tetDonorId, "_", 2)) %>%
  select(libid, full_nt_sequence, v_gene, j_gene, junction, chain, 
         tet, donor, visit) %>%
  pivot_wider(id_cols = c(libid, tet, donor, visit), names_from = chain,
              values_from = c(full_nt_sequence, v_gene, j_gene, junction)) %>%
  tidyr::unnest(cols = c(full_nt_sequence_a, full_nt_sequence_b, v_gene_a, v_gene_b,
                         j_gene_a, j_gene_b, junction_a, junction_b)) %>% 
  group_by(tet, donor, full_nt_sequence_a, full_nt_sequence_b, v_gene_a, v_gene_b,
           j_gene_a, j_gene_b, junction_a, junction_b) %>%
  summarise(count = n()) %>% 
  ungroup() %>%
  filter(!grepl("c\\(", v_gene_b)) %>%
  filter(!is.na(v_gene_b)) %>%
  filter(!grepl("c\\(", v_gene_a)) %>%
  filter(!is.na(v_gene_a)) %>%
  rename(ptid = donor, 
             cdr3_a_nt = full_nt_sequence_a,
             cdr3_b_nt = full_nt_sequence_b,
             v_a_gene = v_gene_a,
             v_b_gene = v_gene_b,
             j_a_gene = j_gene_a,
             j_b_gene = j_gene_b,
             cdr3_a_aa = junction_a,
             cdr3_b_aa = junction_b)

## find triplicates (to compare double alphas)
wen_triplicates <- wen %>%
   group_by(libid) %>%
   filter(n() == 3) %>%
   ungroup()

## drop libids with two betas and make pairs with the double alphas 
double_alphas <- wen_triplicates %>%
   filter(chain == "a" | chain == "b") %>% 
   mutate(tet = str_split_i(tetDonorId, "_", 2)) %>%
   group_by(libid) %>%
   filter(sum(chain == "a") == 2 & sum(chain == "b") < 2) %>%
   select(libid, full_nt_sequence, v_gene, j_gene, junction, chain, 
          tet, donor, visit)

a_chain <- double_alphas %>%
   filter(chain == "a") %>%
   rename_with(~ paste0(., "_a"), -c(libid, tet, donor, visit)) %>%
   select(-chain_a)

b_chain <- double_alphas %>% 
   filter(chain == "b") %>%
   rename_with(~ paste0(., "_b"), -c(libid, tet, donor, visit)) %>%
   select(-chain_b)

paired_data <- full_join(a_chain, b_chain, by = c("libid", "tet", "donor", "visit")) %>%
   group_by(tet, donor, full_nt_sequence_a, full_nt_sequence_b, v_gene_a, v_gene_b,
            j_gene_a, j_gene_b, junction_a, junction_b) %>%
   summarise(count = n()) %>% 
   filter(!is.na(v_gene_b)) %>%
   filter(!is.na(v_gene_a)) %>%
   rename(ptid = donor, 
          cdr3_a_nt = full_nt_sequence_a,
          cdr3_b_nt = full_nt_sequence_b,
          v_a_gene = v_gene_a,
          v_b_gene = v_gene_b,
          j_a_gene = j_gene_a,
          j_b_gene = j_gene_b,
          cdr3_a_aa = junction_a,
          cdr3_b_aa = junction_b) %>%
   mutate(recomb = "TRUE")

pairs$recomb = "FALSE"
all <- rbind (pairs, paired_data)

all$v_b_gene <- paste0(all$v_b_gene, "*01") 
all$j_b_gene <- paste0(all$j_b_gene, "*01")
all$v_a_gene <- paste0(all$v_a_gene, "*01") 
all$j_a_gene <- paste0(all$j_a_gene, "*01")
all$v_a_gene <- gsub("DV", "/DV", all$v_a_gene)

## consolidate by exact nt+aa sequences
check <- all %>%
    group_by(tet, ptid, cdr3_a_nt, cdr3_b_nt, v_a_gene, v_b_gene,
             j_a_gene, j_b_gene, cdr3_a_aa, cdr3_b_aa) %>%
    summarise(count = sum(count))

## consolidate to just unique amino acid sequence
check_aa <- check %>%
  group_by(tet, ptid, cdr3_a_aa, cdr3_b_aa, v_a_gene, v_b_gene,
           j_a_gene, j_b_gene) %>%
  summarise(count = sum(count))
write.csv(check_aa, filename_out, row.names = F)
