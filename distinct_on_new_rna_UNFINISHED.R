## 19/8/22 HAVE NOT FINISHED THIS

#trying out distinct: a method for differential 
#analyses via hierarchical permutation tests

library(tidyverse)
library(distinct)
library(SummarizedExperiment)

source(here::here("create deseq tables from filtered grandslam data.R"))

new_rna_counts_fp <- file.path(here::here("data",'grandslam_new_rna.csv'))
new_rna_counts <- read_csv(new_rna_counts_fp)

grandslam_new_rna_filtered <- new_rna_counts %>% 
  filter(estimate_new_credible_interval < 0.2)

time_1_counts <- create_table_1(grandslam_new_rna_filtered, time = '1')
time_1_meta <- create_metatable(time_1_counts)
time_1_meta = time_1_meta[-1,]
time_1_meta <- rownames_to_column(time_1_meta, var = 'sample_id')

SummarizedExperiment(assays = SimpleList(counts = time_1_counts),
                     colData = time_1_counts,
                     metadata = time_1_meta, 
                     checkDimnames = TRUE)

samples <- time_1_meta$sample_id
group <- time_1_meta$condition
design <- model.matrix(~group)
rownames(design) = samples

res <- distinct_test(x = time_1_counts,
           name_assays_expression = 'logcounts',
           name_cluster = 'condition',
           name_sample = 'sample_id',
           design = design,
           column_to_test = 2,
           min_non_zero_cells = 20,
           n_cores = 2)
