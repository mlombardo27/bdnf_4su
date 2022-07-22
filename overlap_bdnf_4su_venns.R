#4su vs bdnf genes

library(org.Hs.eg.db)
library(tidyverse)
one_v_six <- read_csv('one_vs_6_dds_control.csv')
one_v_six = one_v_six %>% 
  mutate(gene = gsub("\\..*", "", Geneid))

new_ratio_bayesian_p_de <- read_csv("data/new_ratio_bayesian_p_de.csv")
new_ratio_bayesian_p_de = new_ratio_bayesian_p_de %>% 
  dplyr::select(-gene_name) %>% 
  mutate(gene = gsub("\\..*", "", gene)) %>%  
  left_join(annotables::grch38 %>% dplyr::select(ensgene,symbol), by = c("gene" = "ensgene")) %>% 
  unique()

down_2hr_bdnf <- new_ratio_bayesian_p_de %>% 
  filter(time == 2 & total_rna_sig == 'downregulated')

down_4su <- one_v_six %>% 
  filter(log2FoldChange <0 & padj < 0.1)

bdnf_minus_4su_down <- down_2hr_bdnf %>% 
  filter(!(gene %in% down_4su)) %>% 
  write_csv(file = 'just_bdnf_downreg_2hr.csv')

down_4su <- down_4su$gene
down_bdnf <- down_2hr_bdnf$gene

my_list <- list(down_4su = down_4su,
                down_bdnf = down_bdnf)

library(eulerr)

plot(euler(my_list), quantities = TRUE)

#now for upregulated

up_2hr_bdnf <- new_ratio_bayesian_p_de %>% 
  filter(time == 2 & total_rna_sig == 'upregulated')

up_4su <- one_v_six %>% 
  filter(log2FoldChange > 0 & padj < 0.1)

up_4su <- up_4su$gene
up_bdnf <- up_2hr_bdnf$gene

my_list <- list(up_4su = up_4su,
                up_bdnf = up_bdnf)
plot(euler(my_list), quantities = TRUE)

up_bdnf_minus_4su <- up_2hr_bdnf %>% 
  filter(!(gene %in% up_4su)) %>% 
  write_csv(file = 'up_bdnf_minus_4su.csv')
