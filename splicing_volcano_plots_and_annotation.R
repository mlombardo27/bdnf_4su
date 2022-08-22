#splicing volcanos
#code from AL

library(tidyverse)
library(data.table)
library(ggrepel)
source(here::here('splicing_volcanoplot.R'))
# torin_translation = readxl::read_excel("/Users/annaleigh/Downloads/1-s2.0-S1097276519308391-mmc5.xlsx") %>% 
#     janitor::clean_names() %>% 
#     separate_rows(gene_symbol) %>% 
#     as.data.table() 

if (!require("pacman")) install.packages("pacman")
library(pacman)
pacman::p_load(annotatr,janitor,tidyverse,AnnotationDbi,org.Hs.eg.db,plyranges)

#build the annotations for the hg38 basic genes set using the annotatr package
annotations = build_annotations(genome = "hg38", annotations = c("hg38_basicgenes"))
annotations = keepStandardChromosomes(annotations, pruning.mode = "coarse")

# one hour ----------------------------------------------------------------

csv_file_one = fread("/Users/annaleigh/Desktop/bdnf-4su_splicing/control1-bdnf1_annotated_junctions.csv")
one_vp = splicing_dots_tables_function(csv_file_one,list_gene = c("NTRK2","CNOT4","AKT1",'RO60'))
one_vp + ggtitle("Splicing at 1 hour BDNF")

significant_changes_one = csv_file_one %>% 
    filter(probability_changing > 0.9) %>% #90% certain this change is real
    filter(abs(mean_dpsi_per_lsv_junction) > 0.10)  %>%  #Only looking at the events that increase in BNDF 
    dplyr::select(gene_name,gene_id,paste_into_igv_junction, mean_dpsi_per_lsv_junction,
           control1_mean_psi,bdnf1_mean_psi,paste_into_igv_junction,strand,
           seqnames,start,end) %>% 
    mutate(start = start + 1) %>% 
    arrange(-mean_dpsi_per_lsv_junction) %>% 
    makeGRangesFromDataFrame(,keep.extra.columns = TRUE)

my_peak_annotated_one = annotate_regions(regions = significant_changes_one,
                                         annotations = annotations)
my_peak_annotated_one = clean_names(as_tibble(my_peak_annotated_one))

my_peak_annotated_one = my_peak_annotated_one %>% 
    dplyr::select(seqnames,start,end,gene_name,annot_symbol,annot_type,mean_dpsi_per_lsv_junction,
                  control1_mean_psi,bdnf1_mean_psi) %>% 
    unique() %>% 
    filter(!(annot_type %in% c("hg38_genes_1to5kb","hg38_genes_promoters"))) 
    
three_one = my_peak_annotated_one %>% 
    dplyr::select(seqnames,start,end,name,annot_symbol,annot_type,score) %>% 
    unique() %>% 
    filter(!(annot_type %in% c("hg38_genes_1to5kb","hg38_genes_promoters"))) %>% 
    filter(annot_type == "hg38_genes_3UTRs") %>% 
    pull(annot_symbol)
# two hours ---------------------------------------------------------------
csv_file_two = read_csv("C:/Users/mlomb/Desktop/tracked_files_github/bdnf_4su/data/bdnf-4su_splicing/bdnf-4su_splicing/control2-bdnf2_annotated_junctions.csv")
two_vp = splicing_dots_tables_function(csv_file_two,list_gene = c("NTRK2","CNOT4", 'CAMK1', 'SORBS2', 'ABLIM2','PTPN3'))
two_vp + ggtitle("Splicing at 2 hours BDNF")

significant_changes_two = csv_file_two %>% 
    filter(probability_changing > 0.9) %>% #90% certain this change is real
    filter(abs(mean_dpsi_per_lsv_junction) > 0.10)   #Only looking at the events that increase in BNDF 
    dplyr::select(gene_name,gene_id,paste_into_igv_junction, mean_dpsi_per_lsv_junction,
                  bdnf2_mean_psi,control2_mean_psi,paste_into_igv_junction,strand,
                  seqnames,start,end) %>% 
    mutate(start = start + 1) %>% 
    arrange(-mean_dpsi_per_lsv_junction) %>% 
    makeGRangesFromDataFrame(,keep.extra.columns = TRUE)

my_peak_annotated_two = annotate_regions(regions = significant_changes_two,
                                             annotations = annotations)
my_peak_annotated_two = clean_names(as_tibble(my_peak_annotated_two))

my_peak_annotated_two = my_peak_annotated_two %>% 
    dplyr::select(seqnames,start,end,gene_name,annot_symbol,annot_type,
                  mean_dpsi_per_lsv_junction,
                  control2_mean_psi,bdnf2_mean_psi) %>% 
    unique() %>% 
    filter(!(annot_type %in% c("hg38_genes_1to5kb","hg38_genes_promoters"))) 

# six hours ---------------------------------------------------------------
csv_file_six = fread("/Users/annaleigh/Desktop/bdnf-4su_splicing/control6-bdnf6_annotated_junctions.csv")
six_vp = splicing_dots_tables_function(csv_file_six,list_gene = c("NTRK2","CNOT4","AKT1",'RO60'))
six_vp + ggtitle("Splicing at 6 hours BDNF")

significant_changes_six = csv_file_six %>% 
    filter(probability_changing > 0.9) %>% #90% certain this change is real
    filter(abs(mean_dpsi_per_lsv_junction) > 0.10)  %>%  #Only looking at the events that increase in BNDF 
    dplyr::select(gene_name,gene_id,paste_into_igv_junction, mean_dpsi_per_lsv_junction,
                  bdnf6_mean_psi,control6_mean_psi,paste_into_igv_junction,strand,
                  seqnames,start,end) %>% 
    mutate(start = start + 1) %>% 
    arrange(-mean_dpsi_per_lsv_junction) %>% 
    makeGRangesFromDataFrame(,keep.extra.columns = TRUE)

my_peak_annotated_sig_six = annotate_regions(regions = significant_changes_six,
                                             annotations = annotations)
my_peak_annotated_sig_six = clean_names(as_tibble(my_peak_annotated_sig_six))

my_peak_annotated_sig_six = my_peak_annotated_sig_six %>% 
    dplyr::select(seqnames,start,end,gene_name,annot_symbol,
                  annot_type,mean_dpsi_per_lsv_junction,
                  control6_mean_psi,bdnf6_mean_psi) %>% 
    unique() %>% 
    filter(!(annot_type %in% c("hg38_genes_1to5kb","hg38_genes_promoters"))) 


# overlap in time splicing ------------------------------------------------

df_sig_one = csv_file_one %>% 
    filter(probability_changing > 0.9) %>% #90% certain this change is real
    filter(abs(mean_dpsi_per_lsv_junction) > 0.10)  %>%  #Only looking at the events that increase in BNDF 
    dplyr::select(gene_name,gene_id,paste_into_igv_junction, mean_dpsi_per_lsv_junction,
                  control1_mean_psi,bdnf1_mean_psi,paste_into_igv_junction,strand,
                  seqnames,start,end) 

df_sig_two = csv_file_two %>% 
    filter(probability_changing > 0.9) %>% #90% certain this change is real
    filter(abs(mean_dpsi_per_lsv_junction) > 0.10)  %>%  #Only looking at the events that increase in BNDF 
    dplyr::select(gene_name,gene_id,paste_into_igv_junction, mean_dpsi_per_lsv_junction,
                  bdnf2_mean_psi,control2_mean_psi,paste_into_igv_junction,strand,
                  seqnames,start,end) 

df_sig_six = csv_file_six %>% 
    filter(probability_changing > 0.9) %>% #90% certain this change is real
    filter(abs(mean_dpsi_per_lsv_junction) > 0.10)  %>%  #Only looking at the events that increase in BNDF 
    dplyr::select(gene_name,gene_id,paste_into_igv_junction, mean_dpsi_per_lsv_junction,
                  control6_mean_psi,bdnf6_mean_psi,paste_into_igv_junction,strand,
                  seqnames,start,end) 

test = list(splicing_one = unique(df_sig_one$paste_into_igv_junction),
            splicing_two = unique(df_sig_two$paste_into_igv_junction),
            splicing_six = unique(df_sig_six$paste_into_igv_junction))


ggVennDiagram(test, category.names = c("Sig. Splice 1hr",
                                       "Sig. Splice 2hr",
                                       "Sig. Splice 6hr")) + 
    ggtitle("Splice Events 1, 2, 6 hr overlap") +
    scale_color_manual(values = c("black","black","black"))


splice_table = tibble(events = unique(c(df_sig_one$paste_into_igv_junction,
  df_sig_two$paste_into_igv_junction,
  df_sig_six$paste_into_igv_junction))) %>% 
    mutate(one = events %in% df_sig_one$paste_into_igv_junction) %>% 
    mutate(two = events %in% df_sig_two$paste_into_igv_junction) %>% 
    mutate(six = events %in% df_sig_six$paste_into_igv_junction) %>% 
    dplyr::rename(paste_into_igv_junction = events) %>% 
    left_join(rbind(df_sig_one[,.(paste_into_igv_junction,gene_name)],
                    df_sig_two[,.(paste_into_igv_junction,gene_name)],
                    df_sig_six[,.(paste_into_igv_junction,gene_name)])) %>% 
    as.data.table() %>% 
    unique()

#finding genes that change in splicing probability change > 0.9 and deltaPSI >0.2
#seeing what is interesting in terms of RNA decay, transport, BDNF, etc.

# looking at 2hr first

sig_splice_gene_names_2hr <- significant_changes_two %>% 
  dplyr::select(seqnames, gene_name, probability_changing, mean_dpsi_per_lsv_junction,junc_cat) %>% 
  filter(probability_changing > 0.9 & mean_dpsi_per_lsv_junction >0.2) 

#no dynein, dynactin, kinesins genes

#graphing string data

library(ggpubr)

bp_go_2hr_splice_fp <- file.path('data', 'BP GO splicing events 2hr.tsv')
bp_go_2hr_splice <- janitor::clean_names(readr::read_tsv(bp_go_2hr_splice_fp))

bp_go_2hr_splice %>% 
  mutate(term_description = fct_reorder(term_description,strength)) %>% 
  ggplot(aes(y = strength, fill = false_discovery_rate, x = term_description,size = observed_gene_count))+
  geom_point(pch = 21) +
  scale_fill_gradient(low = 'blue', high = 'red')+
  coord_flip()+ 
  theme_pubr()+
  theme(legend.position = 'right')+
  ggtitle("Biological Process GO for Genes with significant splicing changes at 2hr with BDNF")

mf_go_2hr_splice_fp <- file.path('data', 'MF GO splicing events 2hr.tsv')
mf_go_2hr_splice <- janitor::clean_names(readr::read_tsv(mf_go_2hr_splice_fp))

mf_go_2hr_splice %>% 
  mutate(term_description = fct_reorder(term_description,strength)) %>% 
  ggplot(aes(y = strength, fill = false_discovery_rate, x = term_description,size = observed_gene_count))+
  geom_point(pch = 21) +
  scale_fill_gradient(low = 'blue', high = 'red')+
  coord_flip()+ 
  theme_pubr()+
  theme(legend.position = 'right')+
  ggtitle("Molecular Function GO for Genes with significant splicing changes at 2hr with BDNF")

sig_splice_anno_categories_2hr <- GOfuncR::get_anno_categories(sig_splice_gene_names_2hr$gene_name) %>% 

unique(sig_splice_anno_categories_2hr$gene)

names(sig_splice_anno_categories_2hr)[names(sig_splice_anno_categories_2hr) == 'gene'] <- 'gene_name'

significant_changes_two <- significant_changes_two %>% 
  left_join(sig_splice_anno_categories_2hr, by = 'gene_name') %>% 
  filter(probability_changing > 0.9 & mean_dpsi_per_lsv_junction >0.2) %>% 
  dplyr::select(gene_name, mean_dpsi_per_lsv_junction, probability_changing, junc_cat, go_id, name) 

gene_names_significant_changes_two <- unique(significant_changes_two$gene_name)

significant_changes_two %>% 
  filter(grepl('transport|export|signaling|RNA|transcription|splicing|p-body', name, ignore.case = TRUE)) %>%
  pull(gene_name) %>% 
  unique() %>% 
  write.csv('splicing_genes_2hr.csv')


