#breakdown of 4su genes

library(tidyverse)

new_ratio_fp <- file.path(here::here('data','new_ratio_bayesian_p_de.csv'))
new_ratio_de <- read_csv(new_ratio_fp)

new_ratio_de <- new_ratio_de %>% 
  mutate(ensgene = gsub("\\..*", "",gene))

biotypes <- annotables::grch38

new_ratio_de %>% 
  left_join(biotypes, by = 'ensgene') %>% 
  filter(!is.na(biotype)) %>% 
  filter(time == 2 & new_rna_sig == 'bdnf_higher_new_rna' & total_rna_sig == 'not_significant') %>% 
  group_by(biotype) %>% 
  summarize(count = n())%>%
  mutate(biotype = fct_reorder(biotype,count))  %>%
  ggplot(aes(x = biotype, y = count))+
  geom_col()+
  scale_y_continuous(expand = c(0, 0)) +
  ggpubr::theme_pubr()+
  coord_flip()+
  ggtitle("Biotype of Genes", subtitle = 'BDNF Higher Fraction of New RNA with Not Significant Total RNA Change at 2hr')
  
new_ratio_de %>% 
  left_join(biotypes, by = 'ensgene') %>% 
  filter(!is.na(biotype)) %>% 
  filter(time == 2 & new_rna_sig == 'bdnf_higher_new_rna' & total_rna_sig == 'downregulated') %>% 
  group_by(biotype) %>% 
  summarize(count = n())%>%
  mutate(biotype = fct_reorder(biotype,count))  %>%
  ggplot(aes(x = biotype, y = count))+
  geom_col()+
  scale_y_continuous(expand = c(0, 0)) +
  ggpubr::theme_pubr()+
  coord_flip()+
  ggtitle("Biotype of Genes", subtitle = 'BDNF Higher Fraction of New RNA with Downregulated Total RNA Change at 2hr')
