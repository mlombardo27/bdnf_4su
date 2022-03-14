library(tidyverse)
library(PCAtools)
library(gprofiler2)
library(pcaExplorer)
library(clusterProfiler)
library(RColorBrewer)
library(ggpubr)
library(GOfuncR)

grandslam_new_RNA_file_path <- file.path(here::here('data/grandslam_new_rna.csv'))

grandslam_new_rna <- read_csv(grandslam_new_RNA_file_path)

grandslam_new_rna_filtered <- grandslam_new_rna %>% 
  filter(estimate_new_credible_interval < 0.2)

#create counts table for deseq
create_table_1 <- function(filtered_grandslam_rna, 
                           timepoint = "6"){
  
  filtered_df <- filtered_grandslam_rna %>% 
    filter(time == timepoint) %>% 
    select(gene, min2, SampleID) %>% 
    mutate(min2 = round(min2)) %>% 
    pivot_wider(names_from = SampleID, values_from = min2, values_fill = 0) 
  
  return(filtered_df)
}

#create meta table for deseq
create_metatable <- function(table_1){
  
  meta_table <- tibble(sample = table_1 %>%  
                         colnames()) %>% 
    separate(sample, into = c("condition", "well", "time"), remove = FALSE) %>%
    column_to_rownames('sample') %>% 
    mutate(factor_name = factor(condition, levels = c("control","bdnf")))
  
  return(meta_table)
}

#run deseq
library("DESeq2")

min2_deseq <- function(my_table, time){
  
  table_1 <- create_table_1(my_table,time)
  
  table_1 <- table_1 %>% column_to_rownames('gene')
  
  meta_table <- create_metatable(table_1)
  
  dds_min2 <- DESeqDataSetFromMatrix(countData = table_1,
                                     colData = meta_table,
                                     design = ~ factor_name)
  dds_min2 <- DESeq(dds_min2)
  
  
  
  
  res_min2 <- results(dds_min2) %>% 
    as.data.frame() %>% 
    rownames_to_column('ensgene') %>% 
    mutate(ensgene = gsub("\\..*", "", ensgene)) %>%
    left_join(annotables::grch38 %>% select(ensgene,symbol))
  
  return(res_min2)
}


