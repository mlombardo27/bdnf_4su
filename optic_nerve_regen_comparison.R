#looking at this paper about optic nerve regeneration in mice
#seeing if any of our genes are also within this group

library(tidyverse)
library(gprofiler2)

nerve_regen_genes_fp <- file.path("data","optic_nerve_regen_data.xlsx")
nerve_regen_genes <- janitor::clean_names(readxl::read_xlsx(nerve_regen_genes_fp))

nerve_regen_genes <- nerve_regen_genes[-5]

#first have to convert the genes to human

h_genes <- gorth(query = nerve_regen_genes$gene_name, source_organism = 'mmusculus', target_organism = 'hsapiens')

nerve_regen_genes <- nerve_regen_genes %>% 
  left_join(h_genes, by = c('gene_name' = 'input'))

new_ratio_p_fp <- file.path("data",'new_ratio_bayesian_p_de.csv')
new_ratio <- read_csv(new_ratio_p_fp) %>% 
  mutate(gene = gsub("\\..*", "", gene))

stuff_new_ratio <- new_ratio %>% 
  dplyr::select(padj,pvalue, gene, log2FoldChange, new_rna_sig, total_rna_sig, time)

nerve_regen_genes %>% 
  left_join(stuff_new_ratio, by = c('ortholog_ensg' = 'gene')) %>% 
  filter(!is.na(new_rna_sig)) %>% 
  ggplot(aes(x= new_rna_sig, y = on1, fill = total_rna_sig))+
  geom_boxplot()+
  facet_wrap(~time, scales = "free")+
  labs(y = "z-score of change in transport of newly synthesized protein ", title = "Comparing the protein synthesis in optic nerve regeneration with new fraction rna")

