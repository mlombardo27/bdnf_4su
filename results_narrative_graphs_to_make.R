library(tidyverse)

go_genes_2hr_not_sig_total_higher_new_fp <- file.path('data','go_genes_two_hour_not_sig_total_bdnf_higher_new.csv')
go_genes_2hr_not_sig_total_higher_new <- read_csv(go_genes_2hr_not_sig_total_higher_new_fp)

go_genes_2hr_not_sig_total_higher_new <- go_genes_2hr_not_sig_total_higher_new %>% 
  pull(geneID) %>% 
  unique()
  
new_ratio_bayesian_fp <- file.path('data','new_ratio_bayesian_p_de.csv')
new_ratio_bayesian_p_de <- read_csv(new_ratio_bayesian_fp)

smaller_new_ratio_bayesian_p_de <- new_ratio_bayesian_p_de %>% 
  dplyr::select(n_samp_passing_bdnf,n_samp_passing_control,gene_name,time,new_rna_sig, mean_diff, log2FoldChange)

go_genes_i_want <- smaller_new_ratio_bayesian_p_de %>% 
  filter(gene_name %in% go_genes_2hr_not_sig_total_higher_new) %>% 
  filter(time == 2) %>% 
  slice_min(mean_diff, n = 16) %>% 
  pull(gene_name)

go_genes_2hr_downreg_total_higher_new_fp <- file.path('downregulated_twohours_plot_these_genes.csv')
go_genes_2hr_downreg_total_higher_new <- read_csv(go_genes_2hr_downreg_total_higher_new_fp)

go_genes_2hr_downreg_total_higher_new <- go_genes_2hr_downreg_total_higher_new %>% 
  pull(geneID) %>% 
  unique()

go_genes_i_want2 <- smaller_new_ratio_bayesian_p_de %>% 
  filter(gene_name %in% go_genes_2hr_downreg_total_higher_new)  %>% 
  filter(time == 2) %>% 
  slice_min(mean_diff, n = 16) %>% 
  pull(gene_name)

chx_test_genes <- c('DCTN2','DYNC1I2','HMGB2','MDK')

estimate_list_full_fp <- file.path('data','estimate_list_full.csv')
estimate_list_full <- read_csv(estimate_list_full_fp)

test <-estimate_list_full %>% 
  filter(symbol %in% go_genes_i_want) %>% 
  left_join(smaller_new_ratio_bayesian_p_de, by = c('symbol' = 'gene_name', 'time')) %>% 
  unique() %>% 
  filter(n_samp_passing_bdnf == 3 & n_samp_passing_control == 3) %>% 
  mutate(condition = fct_relevel(condition,'control')) %>% 
  mutate(sig_label = ifelse(new_rna_sig %in% c("unclear","bdnf_equals_control"),"", "*")) %>% 
  mutate(sig_label = ifelse(is.na(new_rna_sig),"", sig_label)) %>% 
  group_by(symbol,time) %>% 
  mutate(max_y = max(map)) %>% 
  ungroup()  %>% 
  mutate(max_y = max_y - 0.005)

test2 <-estimate_list_full %>% 
  filter(symbol %in% chx_test_genes) %>% 
  left_join(smaller_new_ratio_bayesian_p_de, by = c('symbol' = 'gene_name', 'time')) %>% 
  unique() %>% 
  filter(n_samp_passing_bdnf == 3 & n_samp_passing_control == 3) %>% 
  mutate(condition = fct_relevel(condition,'control')) %>% 
  mutate(sig_label = ifelse(new_rna_sig %in% c("unclear","bdnf_equals_control"),"", "*")) %>% 
  mutate(sig_label = ifelse(is.na(new_rna_sig),"", sig_label)) %>% 
  group_by(symbol,time) %>% 
  mutate(max_y = max(map)) %>% 
  ungroup()  %>% 
  mutate(max_y = max_y - 0.005)

test3 <-estimate_list_full %>% 
  filter(symbol %in% go_genes_i_want2) %>%  
  left_join(smaller_new_ratio_bayesian_p_de, by = c('symbol' = 'gene_name', 'time')) %>% 
  unique() %>% 
  filter(n_samp_passing_bdnf == 3 & n_samp_passing_control == 3) %>% 
  mutate(condition = fct_relevel(condition,'control')) %>% 
  mutate(sig_label = ifelse(new_rna_sig %in% c("unclear","bdnf_equals_control"),"", "*")) %>% 
  mutate(sig_label = ifelse(is.na(new_rna_sig),"", sig_label)) %>% 
  group_by(symbol,time) %>% 
  mutate(max_y = max(map)) %>% 
  ungroup()  %>% 
  mutate(max_y = max_y - 0.005)

test3 %>% 
  ggplot(aes(x = as.factor(time), y = map, fill = condition))+
  geom_boxplot() +
  geom_point(size = 1.5,pch = 21,position = position_dodge(width = 0.9)) +
  theme(text = element_text(size = 16)) +
  ylab("Fraction New RNA") +
  xlab("Treatment Time") +
  theme(legend.position = 'right') +
  scale_fill_discrete(name = "",labels = c("BDNF","Control")) +
  ggpubr::theme_pubr() +
  theme(legend.position = 'none',
        strip.text.x = element_text( face = "bold.italic"
        )) + 
  scale_fill_manual(values = c('#FF42FF','#0096FF')) +
  scale_color_manual(values = c('#FF42FF','#0096FF')) +
  facet_wrap(~symbol,scales = 'free_y') +
  geom_text(
    inherit.aes = FALSE, (aes(x = as.character(time), y = max_y,label = sig_label)),
    size = 7)+
  ggtitle('Change in the New Fraction of RNA \nBDNF vs Control 2hr \n downregulated total with higher fraction new RNA with BDNF')

new_ratio_bayesian_p_de %>% 
  filter(gene_name %in% go_genes_i_want2) %>% 
  mutate(sig_change = !(padj > 0.1 | is.na(padj))) %>% 
  ggplot(aes(x = time, y= log2FoldChange, fill = sig_change)) +
  geom_hline(yintercept = 0,size = 1)+
  geom_path(aes(x = time, y = log2FoldChange, group = gene_name))+
  geom_point(size = 3,pch = 21)+
  ggtitle('Change in the Total RNA abundance \nBDNF vs Control 2hr') +
  scale_fill_manual(values = c("plum3","turquoise2")) +
  facet_wrap(~gene_name, scales = "free_y")+
  labs(x = "Treatment Time")+
  ggpubr::theme_pubr()+
  theme(legend.position = 'none')

bdnf_early_response_genes <- c('FOS','EGR1','ARC')

new_ratio_bayesian_p_de %>% 
  filter(gene_name %in% bdnf_early_response_genes) %>% 
  mutate(sig_change = !(padj > 0.1 | is.na(padj))) %>% 
  ggplot(aes(x = time, y= log2FoldChange, fill = sig_change)) +
  geom_hline(yintercept = 0,size = 1)+
  geom_path(aes(x = time, y = log2FoldChange, group = gene_name))+
  geom_point(size = 3,pch = 21)+
  ggtitle('Change in the New Fraction of RNA for early response BDNF genes') +
  scale_fill_manual(values = c("darkolivegreen2", "coral")) +
  facet_wrap(~gene_name, scales = "free_y")+
  labs(x = "Treatment Time")+
  ggpubr::theme_pubr()+
  theme(legend.position = 'none')

 #grandslam data: group by sampleid, summarize sum min2, keep time and condition, ungroup, group by time, condition, > mean summed min2

grandslam_fp <- file.path("data",'grandslam_new_rna.csv')
grandslam_new_rna <- read_csv(grandslam_fp)

grandslam_new_rna %>% 
  group_by(SampleID) %>% 
  summarize(sum_min2 = sum(min2), time, condition) %>% 
  ungroup() %>% 
  unique() %>% 
  ggplot(aes(x = condition, y = sum_min2))+
  geom_boxplot()+
  geom_jitter()+
  facet_wrap(~time)+
  scale_y_continuous(labels = scales::comma)+
  labs(y = 'Counts of Converted Reads')

#min2 = the counts of the converted reads that have a minimum of 2 conversions
  
