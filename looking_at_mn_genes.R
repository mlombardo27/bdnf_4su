#looking at motor neuron marker genes identified by this paper to see if bdnf changed the expression

library(tidyverse)

mn_marker_genes_fp <- file.path('data','motorneuron_markers_genes.txt')
mn_marker_genes <- read_lines(mn_marker_genes_fp)

bayesian_fp <- file.path("data","new_ratio_bayesian_p_de.csv")
new_ratio_p <- read_csv(bayesian_fp)

new_ratio_p = new_ratio_p %>% 
  dplyr::select(-gene_name) %>% 
  mutate(gene = gsub("\\..*", "", gene)) %>%  
  left_join(annotables::grch38 %>% dplyr::select(ensgene,symbol), by = c("gene" = "ensgene")) %>% 
  unique()

bdnf_does_not_change_fp <- file.path('bdnf_does_not_change_these_genes_at_2hrs.csv')
bdnf_no_change <- read_csv(bdnf_does_not_change_fp)
idk <- c('UBC')

new_ratio_p %>% 
  filter(symbol %in% mn_marker_genes) %>%
  mutate(sig_change = !(padj > 0.1 | is.na(padj))) %>%
  ggplot(aes(x = time, y = log2FoldChange, fill = sig_change, color = sig_change))+
  geom_hline(yintercept = 0,size = 1)+
  geom_path(aes(x = time, y = log2FoldChange, group = symbol))+
  geom_point()+
  scale_fill_manual(values = c("plum3","turquoise2")) + 
  facet_wrap(~symbol, scales = 'free_y')+
  ggpubr::theme_pubr()+
  ggtitle('Log2FoldChange of Total RNA for Motor Neuron related genes')

#categorizing new rna

full_de_min2 <- readRDS('full_de_min2.rds')

palette<- c("#A6CEE3", "#1F78B4", "#33A02C")

full_de_min2 %>% 
  filter(gene_name %in% mn_marker_genes) %>% 
  group_by(new_rna_sig, total_rna_sig) %>% 
  summarize(countx = n())  %>% 
  ggplot(aes(x = total_rna_sig, y = countx, fill = new_rna_sig))+
  scale_fill_manual(values = palette)+
  geom_col() +
  coord_flip()+
  ggtitle("Categorizing New vs total RNA -\n for Motor Neuron Marker Genes") +
  labs(y = "Number of Genes", 
       x = "Total RNA Categorization",
       fill = "New RNA Categorization")+
  theme(axis.text.y = element_text(angle = 45), legend.position = "bottom", plot.title = element_text(hjust = 0.5), legend.text = element_text(size = 8), legend.title = element_text(size =10))

source(here::here('splicing_volcanoplot.R'))

library(ggrepel)
library(pacman)
pacman::p_load(annotatr,janitor,tidyverse,AnnotationDbi,org.Hs.eg.db,plyranges)

#build the annotations for the hg38 basic genes set using the annotatr package
annotations = build_annotations(genome = "hg38", annotations = c("hg38_basicgenes"))
annotations = keepStandardChromosomes(annotations, pruning.mode = "coarse")  

csv_file_six = read_csv("C:/Users/mlomb/Desktop/tracked_files_github/bdnf_4su/data/bdnf-4su_splicing/bdnf-4su_splicing/control6-bdnf6_annotated_junctions.csv")
six_vp = splicing_dots_tables_function(csv_file_six,list_gene = mn_marker_genes)
six_vp + ggtitle("Splicing at 6 hours BDNF")

csv_file_six %>% 
  filter(gene_name == 'TUBA1A') %>% view()

#finding genes that have no change with bdnf to normalize the qpcr to

new_ratio_p %>% 
  filter(abs(log2FoldChange) < 0.01 & baseMean > 150) %>% 
  filter(time == 2) %>% 
  arrange(-baseMean) %>% 
  write_csv(file = 'bdnf_does_not_change_these_genes_at_2hrs.csv')
