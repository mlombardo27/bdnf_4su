
#results_table this function uses a deseq object's results table,
#the direction will return either the "up" regulated genes or
#the "down" regulated genes.
#this function returns a vector of either up or down regulated
#genes

filter_significant = function(result_table,
                              padj_cutoff = 0.1,
                              log2Fold_cutoff = 0,
                              direction = "up"){
  if(direction == "up"){
    significant_genes <- result_table %>% 
    dplyr::filter(padj < padj_cutoff & log2FoldChange > log2Fold_cutoff) %>% 
      mutate(ensembl = gsub("\\..*", "", Geneid)) %>% 
      arrange(-log2FoldChange) %>% 
      pull(ensembl)
  }
  if(direction == "down"){
    significant_genes <- result_table %>% 
    dplyr::filter(padj < padj_cutoff & log2FoldChange < log2Fold_cutoff) %>% 
      mutate(ensembl = gsub("\\..*", "", Geneid)) %>% 
      arrange(log2FoldChange) %>% 
      pull(ensembl)
  }
  return(significant_genes)
  }
