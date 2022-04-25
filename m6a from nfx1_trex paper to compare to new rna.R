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
  ggpubr::stat_compare_means(label.y = c(2.8,3,2.8))+
  facet_wrap(~new_rna_sig)+
  ggtitle("m6a Motifs on mRNA stability")

s1_table_smol %>% 
  filter(!is.na(time)) %>% 
  unique() %>% 
  ggplot(aes(x = as.character(time), y = m6A_data_sb, fill = total_rna_sig)) +
  geom_boxplot()+
  scale_y_log10()+
  ggpubr::stat_compare_means(label.y = c(2.6,3,2.6))+
  facet_wrap(~new_rna_sig)+
  ggtitle('Expermimentally determined m6a sites on mRNA stability')

s1_table_smol %>% 
  filter(!is.na(time)) %>% 
  ggplot(aes(x = as.character(time), y = GC, fill = total_rna_sig)) +
  geom_boxplot()+
  ggpubr::stat_compare_means(label.y = c(0.8,0.9,0.8))+
  facet_wrap(~new_rna_sig)+
  ggtitle("GC content on mRNA stability")

s1_table_smol %>% 
  filter(!is.na(time)) %>% 
  ggplot(aes(x = as.character(time), y = exonsCount, fill = total_rna_sig)) +
  geom_boxplot()+
  scale_y_log10()+
  ggpubr::stat_compare_means(label.y = c(2.6,2.8,2.6))+
  facet_wrap(~new_rna_sig)+
  ggtitle("Exon Count on mRNA stability")

s1_table_smol %>% 
  filter(!is.na(time)) %>% 
  ggplot(aes(x = as.character(time), y = exonicLength, fill = total_rna_sig)) +
  geom_boxplot()+
  scale_y_log10()+
  ggpubr::stat_compare_means(label.y = c(5,5.5,5))+
  facet_wrap(~new_rna_sig)+
  ggtitle("Exon Length on mRNA stability")

yth_paper_xl_fp <- file.path('data','yth paper excel.xlsx')
yth_paper_xl <- janitor::clean_names(readxl::read_excel(yth_paper_xl_fp))

yth_genes <-  yth_paper_xl %>% 
  filter(category == c('II','II +III')) %>% 
  pull(gene_symbol) %>% 
  unique()

new_rna_genes_2hr_total_rna_downreg <- new_ratio_p %>% 
  filter(time == 2 & new_rna_sig == 'bdnf_higher_new_rna' & total_rna_sig == 'downregulated') %>% 
  pull(gene_name) %>% 
  unique()

list_higher_new_downreg_total <- list(A = yth_genes,
                                      B = new_rna_genes_2hr_total_rna_downreg)

library("ggVennDiagram")

ggVennDiagram(list_higher_new_downreg_total, category.names = c("YTHDF2 targets", "genes with higher new fraction RNA with BDNF and \ndownregulated total rna at 2hr"))


new_rna_genes_2hr_total_rna_not_sig <- new_ratio_p %>% 
  filter(time == 2 & new_rna_sig == 'bdnf_higher_new_rna' & total_rna_sig == 'not_significant') %>% 
  pull(gene_name) %>% 
  unique()

list_higher_new_downreg_total <- list(A = yth_genes,
                                      B = new_rna_genes_2hr_total_rna_not_sig)

ggVennDiagram(list_higher_new_downreg_total, category.names = c("YTHDF2 targets", "genes with higher new fraction RNA with BDNF and \nno total rna change at 2hr"))

new_rna_genes_1hr_total_rna_downreg <- new_ratio_p %>% 
  filter(time == 1 & new_rna_sig == 'bdnf_equals_control' & total_rna_sig == 'downregulated') %>% 
  pull(gene_name) %>% 
  unique()

list_nochange_new_downreg_total <- list(A = yth_genes,
                                      B = new_rna_genes_1hr_total_rna_downreg)

ggVennDiagram(list_nochange_new_downreg_total, category.names = c("YTHDF2 targets", "genes with no change in new rna and \ndownregulated total rna at 1hr"))

new_rna_genes_1hr_total_rna_upreg <- new_ratio_p %>% 
  filter(time == 1 & new_rna_sig == 'bdnf_higher_new_rna' & total_rna_sig == 'upregulated') %>% 
  pull(gene_name) %>% 
  unique()

list_higher_new_upreg_total <- list(A = yth_genes,
                                        B = new_rna_genes_1hr_total_rna_upreg)

ggVennDiagram(list_higher_new_upreg_total, category.names = c("YTHDF2 targets", "genes with higher new rna and \nupregulated total rna at 1hr"))
