res <- run_standard_deseq("C:/Users/mlomb/OneDrive/Desktop/MRes project/Differential Expression/cortical_4su",
                   base_grep = "0",
                   contrast_grep = "24",
                   grep_pattern = "^Ctrl_\\d_(0|24)",
                   baseName = "zero",
                   contrastName = '4suLabelled')

make_volcano_plot(res$results_table)

results <- res$results_table

filtered_res <- filter(results, padj < 0.1 )

ordered_res <- arrange(filtered_res, desc(abs(log2FoldChange)))

top_30 <- head(ordered_res, 30) %>% 
  select(gene_name)

add_name_to_plot(res$results_table, gene_names = top_30$gene_name)

#Now I am doing a GO analysis on this data

library(clusterProfiler)
library(org.Hs.eg.db)
gene <- filtered_res %>% 
  filter(abs(log2FoldChange) > 2) %>% 
  separate(Geneid, into = "ENSEMBL") %>% 
  pull(ENSEMBL)

head(gene)

ggo <- groupGO(gene = gene,
               OrgDb = org.Hs.eg.db,
               ont = "CC",
               keyType = "ENSEMBL",
               level = 3,
               readable = TRUE)

head (ggo)

ego <- enrichGO(gene = gene,
                OrgDb = org.Hs.eg.db,
                ont = "BP",
                keyType = "ENSEMBL",
                pAdjustMethod = "BH",
                pvalueCutoff = 0.01,
                qvalueCutoff = 0.05,
                readable = TRUE)

ggplot(ego)

head(ego)

gene.df <- bitr(gene, fromType = "ENTREZID",
                toType = c("ENSEMBL", "SYMBOL"),
                OrgDb = org.Hs.eg.db)

ego2 <- enrichGO(gene = gene.df$ENSEMBL,
                 OrgDb = org.Hs.eg.db,
                 keyType = 'ENSEMBL',
                 ont = "CC",
                 pAdjustMethod = "BH",
                 pvalueCutoff = 0.01,
                 qvalueCutoff = 0.05)

head(ego2, 3)

ego3 <- gseGO(filtered_res = filtered_res,
              OrgDb = org.Hs.eg.db,
              ont = "CC",
              minGSSize = 100,
              maxGSSize = 500,
              pvalueCutoff = 0.05,
              verbose = FALSE)

goplot(ego)
