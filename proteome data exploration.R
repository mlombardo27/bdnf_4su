#looking at the proteome and the 
#relationships between RNA and protein

library(tidyverse)
library(pcaExplorer)
library(PCAtools)
library(ggpubr)

source(here::here('run_standard_deseq.R'))

total_proteome_fp <- file.path(here::here("data","one_hour_total_proteom.csv"))
total_proteome <- read.csv(total_proteome_fp)

phosphoproteome_fp <- file.path(here::here("data", "one_hour_phosphoproteome.csv"))
phosphoproteome <- read.csv(phosphoproteome_fp)

names(total_proteome)[names(total_proteome) == "target"] <- 'Geneid'

total_rna_feature_counts_fp <- file.path(here::here('data','feature_counts_bdnf_No4'))

one_hour_dds <- run_standard_deseq(total_rna_feature_counts_fp, base_grep = "CONTROL",
                                   contrast_grep = "BDNF", 
                                   grep_pattern = "1hr",
                                   baseName = "control",
                                   contrastName = 'BDNF')

two_hour_dds <- run_standard_deseq(total_rna_feature_counts_fp , base_grep = "CONTROL",
                                   contrast_grep = "BDNF", 
                                   grep_pattern = "2hr",
                                   baseName = "control",
                                   contrastName = 'BDNF')

six_hr_dds <- run_standard_deseq(total_rna_feature_counts_fp , base_grep = "CONTROL",
                                 contrast_grep = "BDNF", 
                                 grep_pattern = "6hr",
                                 baseName = "control",
                                 contrastName = 'BDNF')

time_1_results <- one_hour_dds$results_table
time_2_results <- two_hour_dds$results_table
time_6_results <- six_hr_dds$results_table

de_list <- list(time_1_results,time_2_results,time_6_results)
samp_ids <- c(1,2,6)

full_total_rna_results <-purrr::map2(de_list, samp_ids, ~cbind(.x, time = .y)) %>% 
  data.table::rbindlist(.)

total_rna_results <- full_total_rna_results %>%
  mutate(Geneid = gsub("\\..*", "", Geneid)) %>% 
  left_join(total_proteome, 
            by = c('Geneid'), 
            suffix = c('_total_rna','_proteome'))

new_rna_fp <- file.path(here::here("data",'new_ratio_bayesian_p_de.csv'))
new_rna <- read_csv(new_rna_fp)

 new_rna <- new_rna %>% 
   dplyr::rename(Geneid = gene) %>% 
   mutate(Geneid = gsub("\\..*",'',Geneid))

total_rna_results <- total_rna_results %>% 
  left_join(new_rna,
            by = c('Geneid'),
            suffix = c('_total_rna', '_new_rna')) 

ggplot(total_rna_results, aes(log2FoldChange_total_rna,log_fc, color = new_rna_sig))+
  geom_point(alpha = 0.3)+
  facet_wrap(vars(time_total_rna))

total_rna_results <- total_rna_results %>% 
  mutate(new_rna_sig = case_when(
    new_rna_sig == 'unclear' ~ "Unclear",
    new_rna_sig == 'bdnf_equals_control' ~ "BDNF equals Control",
    new_rna_sig == 'bdnf_higher_new_rna' ~ "BDNF higher new RNA",
    new_rna_sig == 'bdnf_lower_new_rna' ~ "BDNF lower new RNA"
  ))

#total rna pvalue < 0.1
total_rna_results %>% 
  filter(pvalue_total_rna < 0.1) %>% 
  filter(!is.na(new_rna_sig)) %>% 
  ggplot(aes(log2FoldChange_total_rna, log_fc, color = new_rna_sig))+
  geom_point()+
  facet_wrap(~time_total_rna+new_rna_sig)+
  geom_hline(yintercept = 0)+
  geom_vline(xintercept = 0)+
  ggpubr::theme_pubr()+
  labs(color = "New RNA")+
  ylab(expression(paste("Lo", g[2], "Fold Change Proteome")))+
  xlab(expression(paste("Lo", g[2], "Fold Change Total RNA")))+
  ggtitle("Total Proteome (1hr) vs total RNA", subtitle = "Total RNA p-value < 0.1")

#proteome pvalue < 0.1
total_rna_results %>% 
  filter(p_value < 0.1) %>%
  filter(!is.na(new_rna_sig)) %>% 
  ggplot(aes(log2FoldChange_total_rna, log_fc, color = new_rna_sig))+
  geom_point()+
  geom_jitter()+
  facet_wrap(~time_total_rna+new_rna_sig)+
  geom_hline(yintercept = 0)+
  geom_vline(xintercept = 0)+
  ggpubr::theme_pubr()+
  labs(color = "New RNA")+
  ylab(expression(paste("Lo", g[2], "Fold Change Proteome")))+
  xlab(expression(paste("Lo", g[2], "Fold Change Total RNA")))+
  ggtitle("Total Proteome (1hr) vs total RNA", subtitle = "Proteome p-value < 0.1")

#both pvalue < 0.1
total_rna_results %>% 
  filter(p_value < 0.1 & pvalue_total_rna < 0.1) %>% 
  filter(!is.na(new_rna_sig)) %>% 
  ggplot(aes(log2FoldChange_total_rna, log_fc, color = new_rna_sig))+
  geom_point()+
  facet_wrap(~time_total_rna+new_rna_sig)+
  geom_hline(yintercept = 0)+
  geom_vline(xintercept = 0)+
  ggpubr::theme_pubr()+
  labs(color = "New RNA")+
  ylab(expression(paste("Lo", g[2], "Fold Change Proteome")))+
  xlab(expression(paste("Lo", g[2], "Fold Change Total RNA")))+
  ggtitle("Total Proteome (1hr) vs total RNA", subtitle = "Total RNA and Proteome p-values < 0.1")

total_rna_results %>%  
  filter(time_total_rna == 2) %>% 
  filter(log2FoldChange_total_rna > 1) %>% 
  filter(new_rna_sig == 'bdnf_higher_new_rna') %>% 
  filter(log_fc < 0) %>% 
  filter(gene_name_total_rna %in% c('HOMER1', 'SERPINB8', 'PIK3R2', 'GPR50'))
