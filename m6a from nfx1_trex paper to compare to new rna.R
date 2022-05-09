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

#1 higher new, downreg total
new_rna_genes_2hr_total_rna_downreg <- new_ratio_p %>% 
  filter(time == 2 & new_rna_sig == 'bdnf_higher_new_rna' & total_rna_sig == 'downregulated') %>% 
  pull(gene_name) %>% 
  unique()

list_higher_new_downreg_total_2hr <- list(A = yth_genes,
                                      B = new_rna_genes_2hr_total_rna_downreg)

library("ggVennDiagram")

ggVennDiagram(list_higher_new_downreg_total_2hr, category.names = c("YTHDF2 targets", "genes with higher new fraction RNA with BDNF and \ndownregulated total rna at 2hr"))

#2 higher new, not significant total
new_rna_genes_2hr_total_rna_not_sig <- new_ratio_p %>% 
  filter(time == 2 & new_rna_sig == 'bdnf_higher_new_rna' & total_rna_sig == 'not_significant') %>% 
  pull(gene_name) %>% 
  unique()

list_higher_new_not_sig_total_2hr <- list(A = yth_genes,
                                      B = new_rna_genes_2hr_total_rna_not_sig)

ggVennDiagram(list_higher_new_not_sig_total_2hr, category.names = c("YTHDF2 targets", "genes with higher new fraction RNA with BDNF and \nno total rna change at 2hr"))

#3 higher new, total upregulated
new_rna_genes_2hr_total_rna_upreg <- new_ratio_p %>% 
  filter(time == 2 & new_rna_sig == 'bdnf_higher_new_rna' & total_rna_sig == 'upregulated') %>% 
  pull(gene_name) %>% 
  unique()

list_higher_new_upreg_total_2hr <- list(A = yth_genes,
                                      B = new_rna_genes_2hr_total_rna_upreg)

ggVennDiagram(list_higher_new_upreg_total_2hr, category.names = c("YTHDF2 targets", "genes with higher new fraction RNA with BDNF and \nupregulated total rna at 2hr"))

#4 lower new, upregulated total
low_new_genes_2hr_total_rna_upreg <- new_ratio_p %>% 
  filter(time == 2 & new_rna_sig == 'bdnf_lower_new_rna' & total_rna_sig == 'upregulated') %>% 
  pull(gene_name) %>% 
  unique()

list_lower_new_upreg_total_2hr <- list(A = yth_genes,
                                      B = low_new_genes_2hr_total_rna_upreg)

ggVennDiagram(list_lower_new_upreg_total_2hr, category.names = c("YTHDF2 targets", "genes with lower new fraction RNA with BDNF and \nupregulated total rna at 2hr"))

#5 lower new, downregulated total
low_new_genes_2hr_total_rna_downreg <- new_ratio_p %>% 
  filter(time == 2 & new_rna_sig == 'bdnf_lower_new_rna' & total_rna_sig == 'downregulated') %>% 
  pull(gene_name) %>% 
  unique()

list_lower_new_downreg_total_2hr <- list(A = yth_genes,
                                   B = low_new_genes_2hr_total_rna_downreg)

ggVennDiagram(list_lower_new_downreg_total_2hr, category.names = c("YTHDF2 targets", "genes with lower new fraction RNA with BDNF and \ndownregulated total rna at 2hr"))

#6 lower new, not significant total
low_new_genes_2hr_total_rna_not_sig <- new_ratio_p %>% 
  filter(time == 2 & new_rna_sig == 'bdnf_lower_new_rna' & total_rna_sig == 'not_significant') %>% 
  pull(gene_name) %>% 
  unique()

list_lower_new_not_sig_total_2hr <- list(A = yth_genes,
                                     B = low_new_genes_2hr_total_rna_not_sig)

ggVennDiagram(list_lower_new_not_sig_total_2hr, category.names = c("YTHDF2 targets", "genes with lower new fraction RNA with BDNF and \nno sig change in total rna at 2hr"))

#7 not sig new, total not sig
not_sig_new_genes_2hr_total_rna_not_sig <- new_ratio_p %>% 
  filter(time == 2 & new_rna_sig == 'bdnf_equals_control' & total_rna_sig == 'not_significant') %>% 
  pull(gene_name) %>% 
  unique()

list_not_sig_new_not_sig_total_2hr <- list(A = yth_genes,
                                     B = not_sig_new_genes_2hr_total_rna_not_sig)

ggVennDiagram(list_not_sig_new_not_sig_total_2hr, category.names = c("YTHDF2 targets", "genes with not significant new fraction RNA with BDNF and \nno sig change in total rna at 2hr"))

#8 not sig new, upregulated total
not_sig_new_genes_2hr_total_rna_upreg <- new_ratio_p %>% 
  filter(time == 2 & new_rna_sig == 'bdnf_equals_control' & total_rna_sig == 'upregulated') %>% 
  pull(gene_name) %>% 
  unique()

list_not_sig_new_upreg_total_2hr <- list(A = yth_genes,
                                       B = not_sig_new_genes_2hr_total_rna_upreg)

ggVennDiagram(list_not_sig_new_upreg_total_2hr, category.names = c("YTHDF2 targets", "genes with not significant new fraction RNA with BDNF and \nupregulated total rna at 2hr"))

#9 not sig new, downregulated total
not_sig_new_genes_2hr_total_rna_downreg <- new_ratio_p %>% 
  filter(time == 2 & new_rna_sig == 'bdnf_equals_control' & total_rna_sig == 'downregulated') %>% 
  pull(gene_name) %>% 
  unique()

list_not_sig_new_downreg_total_2hr <- list(A = yth_genes,
                                     B = not_sig_new_genes_2hr_total_rna_downreg)

ggVennDiagram(list_not_sig_new_downreg_total_2hr, category.names = c("YTHDF2 targets", "genes with not significant new fraction RNA with BDNF and \ndownregulated total rna at 2hr"))

#1
new_rna_genes_1hr_total_rna_downreg <- new_ratio_p %>% 
  filter(time == 1 & new_rna_sig == 'bdnf_higher_new_rna' & total_rna_sig == 'downregulated') %>% 
  pull(gene_name) %>% 
  unique()

list_higher_new_downreg_total_1hr <- list(A = yth_genes,
                                      B = new_rna_genes_1hr_total_rna_downreg)

ggVennDiagram(list_higher_new_downreg_total_1hr, category.names = c("YTHDF2 targets", "genes with higher new fraction RNA with BDNF and \ndownregulated total rna at 1hr"))


#2
new_rna_genes_1hr_total_rna_not_sig <- new_ratio_p %>% 
  filter(time == 1 & new_rna_sig == 'bdnf_higher_new_rna' & total_rna_sig == 'not_significant') %>% 
  pull(gene_name) %>% 
  unique()

list_higher_new_not_sig_total_1hr <- list(A = yth_genes,
                                      B = new_rna_genes_1hr_total_rna_not_sig)

ggVennDiagram(list_higher_new_not_sig_total_1hr, category.names = c("YTHDF2 targets", "genes with higher new fraction RNA with BDNF and \nno total rna change at 1hr"))

#3
new_rna_genes_1hr_total_rna_upreg <- new_ratio_p %>% 
  filter(time == 1 & new_rna_sig == 'bdnf_higher_new_rna' & total_rna_sig == 'upregulated') %>% 
  pull(gene_name) %>% 
  unique()

list_higher_new_upreg_total_1hr <- list(A = yth_genes,
                                      B = new_rna_genes_1hr_total_rna_upreg)

ggVennDiagram(list_higher_new_upreg_total_1hr, category.names = c("YTHDF2 targets", "genes with higher new fraction RNA with BDNF and \nupregulated total rna at 1hr"))

#4
low_new_genes_1hr_total_rna_upreg <- new_ratio_p %>% 
  filter(time == 1 & new_rna_sig == 'bdnf_lower_new_rna' & total_rna_sig == 'upregulated') %>% 
  pull(gene_name) %>% 
  unique()

list_lower_new_upreg_total_1hr <- list(A = yth_genes,
                                   B = low_new_genes_1hr_total_rna_upreg)

ggVennDiagram(list_lower_new_upreg_total_1hr, category.names = c("YTHDF2 targets", "genes with lower new fraction RNA with BDNF and \nupregulated total rna at 1hr"))

#5
low_new_genes_1hr_total_rna_downreg <- new_ratio_p %>% 
  filter(time == 1 & new_rna_sig == 'bdnf_lower_new_rna' & total_rna_sig == 'downregulated') %>% 
  pull(gene_name) %>% 
  unique()

list_lower_new_downreg_total_1hr <- list(A = yth_genes,
                                     B = low_new_genes_1hr_total_rna_downreg)

ggVennDiagram(list_lower_new_downreg_total_1hr, category.names = c("YTHDF2 targets", "genes with lower new fraction RNA with BDNF and \ndownregulated total rna at 1hr"))

#6
low_new_genes_1hr_total_rna_not_sig <- new_ratio_p %>% 
  filter(time == 1 & new_rna_sig == 'bdnf_lower_new_rna' & total_rna_sig == 'not_significant') %>% 
  pull(gene_name) %>% 
  unique()

list_lower_new_not_sig_total_1hr <- list(A = yth_genes,
                                     B = low_new_genes_1hr_total_rna_not_sig)

ggVennDiagram(list_lower_new_not_sig_total_1hr, category.names = c("YTHDF2 targets", "genes with lower new fraction RNA with BDNF and \nno sig change in total rna at 1hr"))

#7
not_sig_new_genes_1hr_total_rna_not_sig <- new_ratio_p %>% 
  filter(time == 1 & new_rna_sig == 'bdnf_equals_control' & total_rna_sig == 'not_significant') %>% 
  pull(gene_name) %>% 
  unique()

list_not_sig_new_not_sig_total_1hr <- list(A = yth_genes,
                                       B = not_sig_new_genes_1hr_total_rna_not_sig)

ggVennDiagram(list_not_sig_new_not_sig_total_1hr, category.names = c("YTHDF2 targets", "genes with not significant new fraction RNA with BDNF and \nno sig change in total rna at 1hr"))

#8
not_sig_new_genes_1hr_total_rna_upreg <- new_ratio_p %>% 
  filter(time == 1 & new_rna_sig == 'bdnf_equals_control' & total_rna_sig == 'upregulated') %>% 
  pull(gene_name) %>% 
  unique()

list_not_sig_new_upreg_total_1hr <- list(A = yth_genes,
                                     B = not_sig_new_genes_1hr_total_rna_upreg)

ggVennDiagram(list_not_sig_new_upreg_total_1hr, category.names = c("YTHDF2 targets", "genes with not significant new fraction RNA with BDNF and \nupregulated total rna at 1hr"))

#9
not_sig_new_genes_1hr_total_rna_downreg <- new_ratio_p %>% 
  filter(time == 1 & new_rna_sig == 'bdnf_equals_control' & total_rna_sig == 'downregulated') %>% 
  pull(gene_name) %>% 
  unique()

list_not_sig_new_downreg_total_1hr <- list(A = yth_genes,
                                       B = not_sig_new_genes_1hr_total_rna_downreg)

ggVennDiagram(list_not_sig_new_downreg_total_1hr, category.names = c("YTHDF2 targets", "genes with not significant new fraction RNA with BDNF and \ndownregulated total rna at 1hr"))

#1 higher new, downreg total
new_rna_genes_6hr_total_rna_downreg <- new_ratio_p %>% 
  filter(time == 6 & new_rna_sig == 'bdnf_higher_new_rna' & total_rna_sig == 'downregulated') %>% 
  pull(gene_name) %>% 
  unique()

list_higher_new_downreg_total_6hr <- list(A = yth_genes,
                                          B = new_rna_genes_6hr_total_rna_downreg)

ggVennDiagram(list_higher_new_downreg_total_6hr, category.names = c("YTHDF2 targets", "genes with higher new fraction RNA with BDNF and \ndownregulated total rna at 6hr"))

#2 higher new, not significant total
new_rna_genes_6hr_total_rna_not_sig <- new_ratio_p %>% 
  filter(time == 6 & new_rna_sig == 'bdnf_higher_new_rna' & total_rna_sig == 'not_significant') %>% 
  pull(gene_name) %>% 
  unique()

list_higher_new_not_sig_total_6hr <- list(A = yth_genes,
                                          B = new_rna_genes_6hr_total_rna_not_sig)

ggVennDiagram(list_higher_new_not_sig_total_6hr, category.names = c("YTHDF2 targets", "genes with higher new fraction RNA with BDNF and \nno total rna change at 6hr"))

#3 higher new, total upregulated
new_rna_genes_6hr_total_rna_upreg <- new_ratio_p %>% 
  filter(time == 6 & new_rna_sig == 'bdnf_higher_new_rna' & total_rna_sig == 'upregulated') %>% 
  pull(gene_name) %>% 
  unique()

list_higher_new_upreg_total_6hr <- list(A = yth_genes,
                                        B = new_rna_genes_6hr_total_rna_upreg)

ggVennDiagram(list_higher_new_upreg_total_6hr, category.names = c("YTHDF2 targets", "genes with higher new fraction RNA with BDNF and \nupregulated total rna at 6hr"))

#4 lower new, upregulated total
low_new_genes_6hr_total_rna_upreg <- new_ratio_p %>% 
  filter(time == 6 & new_rna_sig == 'bdnf_lower_new_rna' & total_rna_sig == 'upregulated') %>% 
  pull(gene_name) %>% 
  unique()

list_lower_new_upreg_total_6hr <- list(A = yth_genes,
                                       B = low_new_genes_6hr_total_rna_upreg)

ggVennDiagram(list_lower_new_upreg_total_6hr, category.names = c("YTHDF2 targets", "genes with lower new fraction RNA with BDNF and \nupregulated total rna at 6hr"))

#5 lower new, downregulated total
low_new_genes_6hr_total_rna_downreg <- new_ratio_p %>% 
  filter(time == 6 & new_rna_sig == 'bdnf_lower_new_rna' & total_rna_sig == 'downregulated') %>% 
  pull(gene_name) %>% 
  unique()

list_lower_new_downreg_total_6hr <- list(A = yth_genes,
                                         B = low_new_genes_6hr_total_rna_downreg)

ggVennDiagram(list_lower_new_downreg_total_6hr, category.names = c("YTHDF2 targets", "genes with lower new fraction RNA with BDNF and \ndownregulated total rna at 6hr"))

#6 lower new, not significant total
low_new_genes_6hr_total_rna_not_sig <- new_ratio_p %>% 
  filter(time == 6 & new_rna_sig == 'bdnf_lower_new_rna' & total_rna_sig == 'not_significant') %>% 
  pull(gene_name) %>% 
  unique()

list_lower_new_not_sig_total_6hr <- list(A = yth_genes,
                                         B = low_new_genes_6hr_total_rna_not_sig)

ggVennDiagram(list_lower_new_not_sig_total_6hr, category.names = c("YTHDF2 targets", "genes with lower new fraction RNA with BDNF and \nno sig change in total rna at 6hr"))

#7 not sig new, total not sig
not_sig_new_genes_6hr_total_rna_not_sig <- new_ratio_p %>% 
  filter(time == 6 & new_rna_sig == 'bdnf_equals_control' & total_rna_sig == 'not_significant') %>% 
  pull(gene_name) %>% 
  unique()

list_not_sig_new_not_sig_total_6hr <- list(A = yth_genes,
                                           B = not_sig_new_genes_6hr_total_rna_not_sig)

ggVennDiagram(list_not_sig_new_not_sig_total_6hr, category.names = c("YTHDF2 targets", "genes with not significant new fraction RNA with BDNF and \nno sig change in total rna at 6hr"))

#8 not sig new, upregulated total
not_sig_new_genes_6hr_total_rna_upreg <- new_ratio_p %>% 
  filter(time == 6 & new_rna_sig == 'bdnf_equals_control' & total_rna_sig == 'upregulated') %>% 
  pull(gene_name) %>% 
  unique()

list_not_sig_new_upreg_total_6hr <- list(A = yth_genes,
                                         B = not_sig_new_genes_6hr_total_rna_upreg)

ggVennDiagram(list_not_sig_new_upreg_total_6hr, category.names = c("YTHDF2 targets", "genes with not significant new fraction RNA with BDNF and \nupregulated total rna at 6hr"))

#9 not sig new, downregulated total
not_sig_new_genes_6hr_total_rna_downreg <- new_ratio_p %>% 
  filter(time == 6 & new_rna_sig == 'bdnf_equals_control' & total_rna_sig == 'downregulated') %>% 
  pull(gene_name) %>% 
  unique()

list_not_sig_new_downreg_total_6hr <- list(A = yth_genes,
                                           B = not_sig_new_genes_6hr_total_rna_downreg)

ggVennDiagram(list_not_sig_new_downreg_total_6hr, category.names = c("YTHDF2 targets", "genes with not significant new fraction RNA with BDNF and \ndownregulated total rna at 6hr"), show_intersect =TRUE) 

#Doing go analysis on high_new_not_sig_total_2hr, high_new_downreg_total_2hr, high_new_not_sig_total_6hr, high_new_downreg_total_6hr

bp_go_high_new_downreg_total_2hr_fp <- file.path("data", "bp_go_high_new_downreg_total_2hr.tsv")
bp_go_high_new_downreg_total_2hr <- janitor::clean_names(read_tsv(bp_go_high_new_downreg_total_2hr_fp)) %>% 
  mutate(go_type = "biological_process")

mf_go_high_new_downreg_total_2hr_fp <- file.path("data", "mf_go_high_new_downreg_total_2hr.tsv")
mf_go_high_new_downreg_total_2hr <- janitor::clean_names(read_tsv(mf_go_high_new_downreg_total_2hr_fp)) %>% 
  mutate(go_type = 'molecular_function')

cc_go_high_new_downreg_total_2hr_fp <- file.path("data", "cc_go_high_new_downreg_total_2hr.tsv")
cc_go_high_new_downreg_total_2hr <- janitor::clean_names(read_tsv(cc_go_high_new_downreg_total_2hr_fp)) %>% 
  mutate(go_type = 'cellular_component')

combined_go_2hr_high_new_downreg_total <- rbind(bp_go_high_new_downreg_total_2hr,mf_go_high_new_downreg_total_2hr,cc_go_high_new_downreg_total_2hr)

combined_go_2hr_high_new_downreg_total %>%
  slice_max(strength, n = 25) %>%
  mutate(term_description = fct_reorder(term_description,strength)) %>%  
  ggplot(aes(y = strength, fill = false_discovery_rate, x = term_description,size = observed_gene_count)) + 
  geom_point(pch = 21) + 
  facet_wrap(~go_type, scales = 'free')+
  coord_flip()+ 
  theme_pubr()+
  ggtitle("GO of genes overlapping with YTHDF2 from \nhigher new fraction and downregulated total rna at 2hr")

bp_go_high_new_not_sig_total_2hr_fp <- file.path("data", "bp_go_high_new_not_sig_total_2hr.tsv")
bp_go_high_new_not_sig_total_2hr <- janitor::clean_names(read_tsv(bp_go_high_new_not_sig_total_2hr_fp)) %>% 
  mutate(go_type = "biological_process")

mf_go_high_new_not_sig_total_2hr_fp <- file.path("data", "mf_go_high_new_not_sig_total_2hr.tsv")
mf_go_high_new_not_sig_total_2hr <- janitor::clean_names(read_tsv(mf_go_high_new_not_sig_total_2hr_fp)) %>% 
  mutate(go_type = 'molecular_function')

cc_go_high_new_not_sig_total_2hr_fp <- file.path("data", "cc_go_high_new_not_sig_total_2hr.tsv")
cc_go_high_new_not_sig_total_2hr <- janitor::clean_names(read_tsv(cc_go_high_new_not_sig_total_2hr_fp)) %>% 
  mutate(go_type = 'cellular_component')

combined_go_2hr_high_new_not_sig_total <- rbind(bp_go_high_new_not_sig_total_2hr,mf_go_high_new_not_sig_total_2hr,cc_go_high_new_not_sig_total_2hr)

combined_go_2hr_high_new_not_sig_total %>%
  slice_max(strength, n = 20) %>%
  mutate(term_description = fct_reorder(term_description,strength)) %>%  
  ggplot(aes(y = strength, fill = false_discovery_rate, x = term_description,size = observed_gene_count)) + 
  geom_point(pch = 21) + 
  facet_wrap(~go_type, scales = 'free')+
  coord_flip()+ 
  theme_pubr()+
  ggtitle("GO of genes overlapping with YTHDF2 from \nhigher new fraction and not significant total rna at 2hr")

#now 6hr
bp_go_high_new_downreg_total_6hr_fp <- file.path("data", "bp_go_high_new_downreg_total_6hr.tsv")
bp_go_high_new_downreg_total_6hr <- janitor::clean_names(read_tsv(bp_go_high_new_downreg_total_6hr_fp)) %>% 
  mutate(go_type = "biological_process")

mf_go_high_new_downreg_total_6hr_fp <- file.path("data", "mf_go_high_new_downreg_total_6hr.tsv")
mf_go_high_new_downreg_total_6hr <- janitor::clean_names(read_tsv(mf_go_high_new_downreg_total_6hr_fp)) %>% 
  mutate(go_type = 'molecular_function')

cc_go_high_new_downreg_total_6hr_fp <- file.path("data", "cc_go_high_new_downreg_total_6hr.tsv")
cc_go_high_new_downreg_total_6hr <- janitor::clean_names(read_tsv(cc_go_high_new_downreg_total_6hr_fp)) %>% 
  mutate(go_type = 'cellular_component')

combined_go_6hr_high_new_downreg_total <- rbind(bp_go_high_new_downreg_total_6hr,mf_go_high_new_downreg_total_6hr,cc_go_high_new_downreg_total_6hr)

combined_go_6hr_high_new_downreg_total %>%
  slice_max(strength, n = 25) %>%
  mutate(term_description = fct_reorder(term_description,strength)) %>%  
  ggplot(aes(y = strength, fill = false_discovery_rate, x = term_description,size = observed_gene_count)) + 
  geom_point(pch = 21) + 
  facet_wrap(~go_type, scales = 'free')+
  coord_flip()+ 
  theme_pubr()+
  ggtitle("GO of genes overlapping with YTHDF2 from \nhigher new fraction and downregulated total rna at 6hr")

bp_go_high_new_not_sig_total_6hr_fp <- file.path("data", "bp_go_high_new_not_sig_total_6hr.tsv")
bp_go_high_new_not_sig_total_6hr <- janitor::clean_names(read_tsv(bp_go_high_new_not_sig_total_6hr_fp)) %>% 
  mutate(go_type = "biological_process")

mf_go_high_new_not_sig_total_6hr_fp <- file.path("data", "mf_go_high_new_not_sig_total_6hr.tsv")
mf_go_high_new_not_sig_total_6hr <- janitor::clean_names(read_tsv(mf_go_high_new_not_sig_total_6hr_fp)) %>% 
  mutate(go_type = 'molecular_function')

cc_go_high_new_not_sig_total_6hr_fp <- file.path("data", "cc_go_high_new_not_sig_total_6hr.tsv")
cc_go_high_new_not_sig_total_6hr <- janitor::clean_names(read_tsv(cc_go_high_new_not_sig_total_6hr_fp)) %>% 
  mutate(go_type = 'cellular_component')

combined_go_6hr_high_new_not_sig_total <- rbind(bp_go_high_new_not_sig_total_6hr,mf_go_high_new_not_sig_total_6hr,cc_go_high_new_not_sig_total_6hr)

combined_go_6hr_high_new_not_sig_total %>%
  slice_max(strength, n = 20) %>%
  mutate(term_description = fct_reorder(term_description,strength)) %>%  
  ggplot(aes(y = strength, fill = false_discovery_rate, x = term_description,size = observed_gene_count)) + 
  geom_point(pch = 21) + 
  facet_wrap(~go_type, scales = 'free')+
  coord_flip()+ 
  theme_pubr()+
  ggtitle("GO of genes overlapping with YTHDF2 from \nhigher new fraction and not significant total rna at 6hr")
