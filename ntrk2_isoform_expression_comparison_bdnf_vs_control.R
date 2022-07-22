# im taking the salmon data and looking at the expression levels of ntrk2 for each timepoint
#below creates boxplots comparing isoform expression of ntrk2 for bdnf vs control at 1,2,6hrs

file_path <- "C:/Users/mlomb/Desktop/tracked_files_github/bdnf_4su/data/Full BDNF 4SU - salmon/Full BDNF 4SU - salmon/one_hour"
suffix <- "quant.sf"
estimate_files <- list.files(file_path,
                             recursive = TRUE,
                            pattern = 'quant.sf',
                            full.names = TRUE)

salmon_full = purrr::map(estimate_files,read_tsv)

samp_ids <- gsub(estimate_files,pattern = file_path, replacement = '')
samp_ids <- gsub(samp_ids, pattern = '/', replacement = '')
samp_ids <- gsub(samp_ids, pattern = 'quant.sf', replacement = '')
  
##add on the name as an additional column
salmon_full = purrr::map2(salmon_full, samp_ids, ~cbind(.x, sample = .y))
salmon_full = data.table::rbindlist(salmon_full)

salmon_full <- salmon_full %>% 
  mutate(Name = gsub('\\..*',"", Name)) %>% 
  left_join(annotables::grch38_tx2gene, by = c(Name = 'enstxp')) %>% 
  left_join(annotables::grch38 %>%  dplyr::select(symbol, ensgene), by = 'ensgene')

ntrk2 <- salmon_full %>% 
  filter(symbol == 'NTRK2')

ntrk2 <- ntrk2 %>% 
  mutate(condition = ifelse(grepl(x=sample, 'CONTROL'), 'CONTROL', 'BDNF'))

ntrk2 %>% 
  filter(!(Name == 'ENST00000376208')) %>% 
  ggplot(aes(x = Name, y = TPM, fill = condition))+
  geom_boxplot()+
  coord_flip()+
  ggpubr::theme_pubclean()+
  ggpubr::stat_compare_means(label = 'p.signif')+
  ggtitle('Number of reads (in TPM) for each of the NTRK2 isoforms at 1 hour')+
  geom_point(pch = 21,position = position_jitterdodge(jitter.width = 0.3), size = 2) 

ntrk2 %>% 
  ggplot(aes(x = Name, y = NumReads))+
  geom_boxplot()+
  coord_flip()+
  ggpubr::theme_pubclean()+
  ggpubr::stat_compare_means(label = 'p.signif')+
  ggtitle('Number of reads (in raw reads) for each of the NTRK2 isoforms at 1 hour')

#2hr

file_path <- "C:/Users/mlomb/Desktop/tracked_files_github/bdnf_4su/data/Full BDNF 4SU - salmon/Full BDNF 4SU - salmon/two_hour"
suffix <- "quant.sf"
estimate_files <- list.files(file_path,
                             recursive = TRUE,
                             pattern = 'quant.sf',
                             full.names = TRUE)

salmon_full = purrr::map(estimate_files,read_tsv)

samp_ids <- gsub(estimate_files,pattern = file_path, replacement = '')
samp_ids <- gsub(samp_ids, pattern = '/', replacement = '')
samp_ids <- gsub(samp_ids, pattern = 'quant.sf', replacement = '')

##add on the name as an additional column
salmon_full = purrr::map2(salmon_full, samp_ids, ~cbind(.x, sample = .y))
salmon_full = data.table::rbindlist(salmon_full)

salmon_full <- salmon_full %>% 
  mutate(Name = gsub('\\..*',"", Name)) %>% 
  left_join(annotables::grch38_tx2gene, by = c(Name = 'enstxp')) %>% 
  left_join(annotables::grch38 %>%  dplyr::select(symbol, ensgene), by = 'ensgene')

ntrk2 <- salmon_full %>% 
  filter(symbol == 'NTRK2')

ntrk2 <- ntrk2 %>% 
  mutate(condition = ifelse(grepl(x=sample, 'CONTROL'), 'CONTROL', 'BDNF'))

ntrk2 %>% 
  ggplot(aes(x = Name, y = TPM, fill = condition))+
  geom_boxplot()+
  coord_flip()+
  ggpubr::theme_pubclean()+
  ggpubr::stat_compare_means(label = 'p.signif')+
  ggtitle('Number of reads (in TPM) for each of the NTRK2 isoforms at 2 hours')+
  geom_point(pch = 21,position = position_jitterdodge(jitter.width = 0.3), size = 2) 

# 6hr

file_path <- "C:/Users/mlomb/Desktop/tracked_files_github/bdnf_4su/data/Full BDNF 4SU - salmon/Full BDNF 4SU - salmon/six_hour"
suffix <- "quant.sf"
estimate_files <- list.files(file_path,
                             recursive = TRUE,
                             pattern = 'quant.sf',
                             full.names = TRUE)

salmon_full = purrr::map(estimate_files,read_tsv)

samp_ids <- gsub(estimate_files,pattern = file_path, replacement = '')
samp_ids <- gsub(samp_ids, pattern = '/', replacement = '')
samp_ids <- gsub(samp_ids, pattern = 'quant.sf', replacement = '')

##add on the name as an additional column
salmon_full = purrr::map2(salmon_full, samp_ids, ~cbind(.x, sample = .y))
salmon_full = data.table::rbindlist(salmon_full)

salmon_full <- salmon_full %>% 
  mutate(Name = gsub('\\..*',"", Name)) %>% 
  left_join(annotables::grch38_tx2gene, by = c(Name = 'enstxp')) %>% 
  left_join(annotables::grch38 %>%  dplyr::select(symbol, ensgene), by = 'ensgene')

ntrk2 <- salmon_full %>% 
  filter(symbol == 'NTRK2')

ntrk2 <- ntrk2 %>% 
  mutate(condition = ifelse(grepl(x=sample, 'CONTROL'), 'CONTROL', 'BDNF'))

ntrk2 %>% 
  ggplot(aes(x = Name, y = TPM, fill = condition))+
  geom_boxplot()+
  coord_flip()+
  ggpubr::theme_pubclean()+
  ggpubr::stat_compare_means(label = 'p.signif')+
  ggtitle('Number of reads (in TPM) for each of the NTRK2 isoforms at 6 hours')+
  geom_point(pch = 21,position = position_jitterdodge(jitter.width = 0.3),size = 2) 


         
         