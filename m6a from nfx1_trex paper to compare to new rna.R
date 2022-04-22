library(tidyverse)

s1_table_fp <- file.path(here::here('data','table_s1_nfx_trex_paper.xlsx'))
s1_table <- readxl::read_excel(s1_table_fp)

new_ratio_fp <- file.path(here::here('data','new_ratio_bayesian_p_de.csv'))
new_ratio_p <- read_csv(new_ratio_fp)

s1_table_smol <- s1_table %>% 
  select(gene_name, GC, exonsCount, exonicLength, m6A, m6A_data_sb) %>% 
  left_join(new_ratio_p, by = 'gene_name')

s1_table_smol %>% 
  filter(!is.na(time)) %>% 
  ggplot(aes(x = as.character(time), y = m6A, fill = total_rna_sig)) +
  geom_boxplot()+
  ggpubr::stat_compare_means()+
  facet_wrap(~new_rna_sig)

s1_table_smol %>% 
  filter(!is.na(time)) %>% 
  unique() %>% 
  ggplot(aes(x = as.character(time), y = m6A_data_sb, fill = total_rna_sig)) +
  geom_boxplot()+
  ggpubr::stat_compare_means()+
  facet_wrap(~new_rna_sig)
