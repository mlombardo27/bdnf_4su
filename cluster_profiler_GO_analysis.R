#this function will take in the up or down regulated vector from
#filter_significant and return a list of go enrichment for
#the process of interest (bp, cc, mf)

cp_GO_analysis = function(significant_genes, 
                          bg_ids,
                          typeGO = "cc"){
  library(tidyverse)
  library(clusterProfiler)
  
  # bg_ids <- dds_results$results_table$Geneid %>% 
  #   gsub("\\..*", "", .)
  
  if(typeGO == "cc"){
    go_result <- enrichGO(gene          = significant_genes,
             universe      = bg_ids,
             keyType = "ENSEMBL",
             OrgDb         = org.Hs.eg.db,
             ont           = typeGO,
             pAdjustMethod = "BH",
             pvalueCutoff  = 0.01,
             qvalueCutoff  = 0.05,
             readable      = TRUE)
  }
  
  if(typeGO == "bp"){
    go_result <- enrichGO(gene          = significant_genes,
             universe      = bg_ids,
             keyType = "ENSEMBL",
             OrgDb         = org.Hs.eg.db,
             ont           = typeGO,
             pAdjustMethod = "BH",
             pvalueCutoff  = 0.01,
             qvalueCutoff  = 0.05,
             readable      = TRUE)
  }
  if(typeGO == "mf"){
    go_result <- enrichGO(gene          = significant_genes,
             universe      = bg_ids,
             keyType = "ENSEMBL",
             OrgDb         = org.Hs.eg.db,
             ont           = typeGO,
             pAdjustMethod = "BH",
             pvalueCutoff  = 0.01,
             qvalueCutoff  = 0.05,
             readable      = TRUE)
  }
  return(go_result)
}

