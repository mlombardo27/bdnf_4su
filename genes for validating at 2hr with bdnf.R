#Jobert wanted to validate some of the genes so I selected genes that were
#either downregulated/not significant with a high level of new rna from BDNF 
#at 2hrs 
#we're going to validate these using qPCR
#here's the criteria

gene_2hr_high4su_downreg <- new_ratio_p_de %>% 
  filter(time ==2, new_rna_sig == 'bdnf_higher_new_rna', total_rna_sig == 'downregulated') %>% 
  slice_min(log2FoldChange, n=10) %>% 
  arrange(-mean_diff) %>% 
  select(gene, gene_name, mean_diff, log2FoldChange, mean_bdnf_ntr, mean_control_ntr)

write_csv(gene_2hr_high4su_downreg, "top10gene_2hr_high4su_downreg.csv")

gene_2hr_high4su_notsig <- new_ratio_p_de %>% 
  filter(time ==2, new_rna_sig == 'bdnf_higher_new_rna', total_rna_sig == 'not_significant') %>%
  filter(!is.na(log2FoldChange)) %>% 
  slice_min(mean_diff, n=10) %>% 
  arrange(-mean_bdnf_ntr) %>% 
  select(gene, gene_name, mean_diff, log2FoldChange, mean_bdnf_ntr, mean_control_ntr)

write_csv(gene_2hr_high4su_notsig, "top10gene_2hr_notsig.csv")


