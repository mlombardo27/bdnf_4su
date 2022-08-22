library(tidyverse)

up_regulated = function(res){
  res <- res %>% 
    dplyr::filter(padj < 0.05 & log2FoldChange > 0 ) %>% 
  mutate(ensembl = gsub("\\..*", "", Geneid)) %>% 
  arrange(-log2FoldChange) %>% 
  pull(ensembl)
}
  
results_1hr <- one_hour_dds$results_table

up_1hr_ensembl <- results_1hr %>% 
  dplyr::filter(padj < 0.05 & log2FoldChange > 0 ) %>% 
  mutate(ensembl = gsub("\\..*", "", Geneid)) %>% 
  arrange(-log2FoldChange) %>% 
  pull(ensembl)