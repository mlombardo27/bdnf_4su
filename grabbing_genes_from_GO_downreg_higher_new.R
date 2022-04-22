
cc_go_2hr_downreg_high_new <- janitor::clean_names(read_tsv(here::here('data','cc 2 downreg bdnf higher new.tsv')))

bp_go_2hr_downreg_high_new <- janitor::clean_names(read_tsv(here::here('data', 'BP 2 downreg bdnf higher new .tsv')))

mf_go_2hr_downreg_high_new <- janitor::clean_names(read_tsv(here::here('data', 'mf 2 downreg bdnf higher new.tsv')))

genes_of_interest_2hr_bp <-bp_go_2hr_downreg_high_new %>% 
  filter(grepl('splic', term_description)) %>% 
  pull(matching_proteins_in_your_network_labels) 

genes_of_interest_2hr_bp <- genes_of_interest_2hr_bp %>% 
  str_split(',') %>% 
  purrr::simplify() %>% 
  unique()

genes_of_interest_2hr_cc <-cc_go_2hr_downreg_high_new %>% 
  filter(grepl('splic', term_description, ignore.case = TRUE)) %>% 
  pull(matching_proteins_in_your_network_labels) 
         
genes_of_interest_2hr_cc <- genes_of_interest_2hr_cc %>% 
  str_split(',') %>% 
  purrr::simplify() %>% 
  unique()

genes_of_interest_2hr_mf <-mf_go_2hr_downreg_high_new %>% 
  filter(grepl('splic', term_description, ignore.case = TRUE)) %>% 
  pull(matching_proteins_in_your_network_labels) 
#no splicing terms

all_genes_of_interest_2hr <- c(genes_of_interest_2hr_bp, genes_of_interest_2hr_cc) %>% 
  unique() %>% 
  clipr::write_clip()


