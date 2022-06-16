#chisquare on isoform switching

library(tidyverse)
library(IsoformSwitchAnalyzeR)
library(corrplot)

new_ratio_bayesian_fp <- file.path('data','new_ratio_bayesian_p_de.csv')
new_ratio <- read_csv(new_ratio_bayesian_fp)

isoform_switch <- readRDS('isoform_switch_2hr_data.rds')

isoform_switch_consequence <- isoform_switch$switchConsequence

coding_to_noncoding_genes <- read_csv('coding_to_noncoding_genes_2hr.csv')

#on both new and total rna categories
my_table <- new_ratio %>% 
  filter(time == 2) %>% 
  dplyr::select(gene_name, new_rna_sig, total_rna_sig) %>% 
  left_join(isoform_switch_consequence, by = 'gene_name') %>% 
  mutate(switchConsequence = ifelse(is.na(switchConsequence), 'no_switch', switchConsequence)) %>% 
  mutate(new_total_sig_cat = glue::glue("{new_rna_sig}_{total_rna_sig}"))

my_table <- my_table %>% 
  group_by(switchConsequence, new_total_sig_cat) %>% 
  summarize(cnt = n()) %>% 
  ungroup()

my_table <- my_table %>% 
  pivot_wider(names_from = 'new_total_sig_cat', values_from = 'cnt', values_fill = 0) %>% 
  column_to_rownames('switchConsequence')

chisq <- chisq.test(my_table)

corrplot(chisq$residuals, is.cor = FALSE)

contrib <- 100*chisq$residuals^2/chisq$statistic %>% 
  round(3)

corrplot(contrib, is.cor = FALSE)

#on total rna category
my_table <- new_ratio %>% 
  filter(time == 2) %>% 
  dplyr::select(gene_name, new_rna_sig, total_rna_sig) %>% 
  left_join(isoform_switch_consequence, by = 'gene_name') %>% 
  mutate(switchConsequence = ifelse(is.na(switchConsequence), 'no_switch', switchConsequence)) 

my_table <- my_table %>% 
  group_by(switchConsequence, total_rna_sig) %>% 
  summarize(cnt = n()) %>% 
  ungroup()

my_table <- my_table %>% 
  pivot_wider(names_from = 'total_rna_sig', values_from = 'cnt', values_fill = 0) %>% 
  column_to_rownames('switchConsequence')

chisq <- chisq.test(my_table)

corrplot(chisq$residuals, is.cor = FALSE)

contrib <- 100*chisq$residuals^2/chisq$statistic %>% 
  round(3)

corrplot(contrib, is.cor = FALSE)

#on new rna category
my_table <- new_ratio %>% 
  filter(time == 2) %>% 
  dplyr::select(gene_name, new_rna_sig, total_rna_sig) %>% 
  left_join(isoform_switch_consequence, by = 'gene_name') %>% 
  mutate(switchConsequence = ifelse(is.na(switchConsequence), 'no_switch', switchConsequence)) 

my_table <- my_table %>% 
  group_by(switchConsequence, new_rna_sig) %>% 
  summarize(cnt = n()) %>% 
  ungroup()

my_table <- my_table %>% 
  pivot_wider(names_from = 'new_rna_sig', values_from = 'cnt', values_fill = 0) %>% 
  column_to_rownames('switchConsequence')

chisq <- chisq.test(my_table)

corrplot(chisq$residuals, is.cor = FALSE)

contrib <- 100*chisq$residuals^2/chisq$statistic %>% 
  round(3)

corrplot(contrib, is.cor = FALSE)

isoform_switch_consequence %>% 
  filter(switchConsequence == 'Intron retention switch')

switchPlot(isoform_switch, gene = 'METTL23')

