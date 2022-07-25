library(tidyverse)

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

# plotting ntrk2 isoforms from g38
library(magrittr)
library(ggtranscript)

g38 <- readRDS('gencode.v38.gtf.RDS')

ntrk2_from_g38 <- g38 %>% 
  filter(gene_name == 'NTRK2') %>% 
  select(seqnames, start, end, strand, type, gene_name, gene_id, transcript_id, transcript_name, transcript_type, exon_number)

ntrk2_i_care_about <- ntrk2_from_g38 %>% 
  mutate(transcript_id = gsub('\\..*',"", transcript_id)) %>% 
  filter(transcript_id == c('ENST00000277120','ENST00000323115','ENST00000376213','ENST00000359847', 'ENST00000304053'))

ntrk2_cds <- ntrk2_i_care_about %>% 
  filter(type == 'CDS')

ntrk2_exons <- ntrk2_i_care_about %>% 
  filter(type == 'exon')

ntrk2_exons_prot_coding <- ntrk2_exons %>% 
  filter(transcript_type == 'protein_coding')

types_i_want <- ntrk2_i_care_about %>% filter(!(type =='transcript')) %>% pull(type) %>% unique()

ntrk2_i_care_about %>% 
  ggplot(aes(xstart = start,
             xend = end,
             y = transcript_id))+
  geom_range(aes(fill = type))+
  geom_range(data = ntrk2_cds)+
  geom_intron(data = to_intron(ntrk2_i_care_about, 'transcript_name'),
              aes(strand = strand),
              arrow.min.intron.length = 500)

ntrk2_120 <- ntrk2_i_care_about %>% filter(transcript_id == 'ENST00000277120')
ntrk2_115<- ntrk2_i_care_about %>% filter(transcript_id == 'ENST00000323115')

ntrk2_115 %>% 
  ggplot(aes(xstart = start,
             xend = end,
             y = 'NTRK2-115/120',
             fill = type))+
  geom_half_range(aes(fill = type))+
  geom_intron(data = to_intron(ntrk2_115, 'transcript_name'),
              arrow.min.intron.length = 500)+
  geom_half_range(data = ntrk2_120, 
                  range.orientation = 'top')+
  geom_intron(data = to_intron(ntrk2_120, 'transcript_name'),
              arrow.min.intron.length = 500)+
  coord_cartesian(xlim = c(84670000, 84720000))

ntrk2_rescaled <- shorten_gaps(exons = ntrk2_exons,
                              introns = to_intron(ntrk2_exons, 'transcript_name'),
                              group_var = 'transcript_name')

ntrk2_rescaled_exons <- ntrk2_rescaled %>% filter(type == 'exon')
ntrk2_rescaled_introns <- ntrk2_rescaled %>% filter(type == 'intron')

ntrk2_rescaled_exons_prot_cod <- ntrk2_rescaled_exons %>% 
  filter(transcript_type =='protein_coding')

ntrk2_rescaled_exons_prot_cod %>% 
  ggplot(aes(
    xstart = start,
    xend = end,
    y = transcript_name)) +
  geom_range(aes(height = 0.25))+
  geom_intron(
    data = ntrk2_rescaled_introns,
    aes(strand = strand), 
    arrow.min.intron.length = 300)

# going to look at the junctions for ntrk2 isoforms 213, 115, 120 that are getting changed the most with BDNF

splicing_junctions_2hr <- read_csv("C:/Users/mlomb/Desktop/tracked_files_github/bdnf_4su/data/bdnf-4su_splicing/bdnf-4su_splicing/control2-bdnf2_annotated_junctions.csv")

#for transcript 115
ntrk2_115_exons <- ntrk2_exons %>% 
  filter(transcript_id == 'ENST00000323115')
ntrk2_115_cds <- ntrk2_cds %>% 
  filter(transcript_id == 'ENST00000323115')

ntrk2_junctions <- splicing_junctions_2hr %>% 
  filter(gene_name == 'NTRK2') %>% 
  mutate(transcript_id = 'ENST00000323115')

ntrk2_115_exons %>% 
  ggplot(aes(xstart = start,
             xend= end,
             y = transcript_id))+
  geom_range(data = ntrk2_115_cds, aes(fill = type))+
  geom_range(fill = 'white',
             height = 0.25)+
  geom_range(data = ntrk2_115_cds)+
  geom_intron(data = to_intron(ntrk2_115_exons, 'transcript_name'))+
  geom_junction(data = ntrk2_junctions,
                aes(size = bdnf2_mean_psi),
                junction.y.max = 0.5)+
  geom_junction_label_repel(data = ntrk2_junctions,
                            aes(label = round(bdnf2_mean_psi, 2)),
                            junction.y.max = 0.5)+
  scale_size_continuous(range = c(0.1,1))+
  coord_cartesian(xlim = c(84667500, 84712000))

#for transcript 120-- most expressed transcript

ntrk2_120_exons <- ntrk2_exons %>% 
  filter(transcript_id == 'ENST00000277120')
ntrk2_120_cds <- ntrk2_cds %>% 
  filter(transcript_id == 'ENST00000277120')

ntrk2_junctions <- splicing_junctions_2hr %>% 
  filter(gene_name == 'NTRK2') %>% 
  mutate(transcript_id = 'ENST00000277120')

ntrk2_120_exons %>% 
  ggplot(aes(xstart = start,
             xend= end,
             y = transcript_id))+
  geom_range(data = ntrk2_120_cds, aes(fill = type))+
  geom_range(fill = 'white',
             height = 0.25)+
  geom_range(data = ntrk2_120_cds)+
  geom_intron(data = to_intron(ntrk2_120_exons, 'transcript_name'))+
  geom_junction(data = ntrk2_junctions,
                aes(size = bdnf2_mean_psi),
                junction.y.max = 0.5)+
  geom_junction_label_repel(data = ntrk2_junctions,
                            aes(label = round(bdnf2_mean_psi, 2)),
                            junction.y.max = 0.5)+
  scale_size_continuous(range = c(0.1,1))+
  coord_cartesian(xlim = c(84667500, 84675000))

# 213--changes with control vs bdnf

ntrk2_213_exons <- ntrk2_exons %>% 
  filter(transcript_id == 'ENST00000376213')
ntrk2_213_cds <- ntrk2_cds %>% 
  filter(transcript_id == 'ENST00000376213')

ntrk2_junctions <- splicing_junctions_2hr %>% 
  filter(gene_name == 'NTRK2') %>% 
  mutate(transcript_id = 'ENST00000376213')

ntrk2_213_exons %>% 
  ggplot(aes(xstart = start,
             xend= end,
             y = transcript_id))+
  geom_range(data = ntrk2_213_cds, aes(fill = type))+
  geom_range(fill = 'white',
             height = 0.25)+
  geom_range(data = ntrk2_213_cds)+
  geom_intron(data = to_intron(ntrk2_213_exons, 'transcript_name'))+
  geom_junction(data = ntrk2_junctions,
                aes(size = bdnf2_mean_psi),
                junction.y.max = 0.5)+
  geom_junction_label_repel(data = ntrk2_junctions,
                            aes(label = round(bdnf2_mean_psi, 2)),
                            junction.y.max = 0.5)+
  scale_size_continuous(range = c(0.1,1))+
  coord_cartesian(xlim = c(84668000, 84705000))

#doing some extra ones 
# 847- changes with BDNF(signficant)

ntrk2_847_exons <- ntrk2_exons %>% 
  filter(transcript_id == 'ENST00000359847')
ntrk2_847_cds <- ntrk2_cds %>% 
  filter(transcript_id == 'ENST00000359847')

ntrk2_junctions <- splicing_junctions_2hr %>% 
  filter(gene_name == 'NTRK2') %>% 
  mutate(transcript_id = 'ENST00000359847')

ntrk2_847_exons %>% 
  ggplot(aes(xstart = start,
             xend= end,
             y = transcript_id))+
  geom_range(data = ntrk2_847_cds, aes(fill = type))+
  geom_range(fill = 'white',
             height = 0.25)+
  geom_range(data = ntrk2_847_cds)+
  geom_intron(data = to_intron(ntrk2_847_exons, 'transcript_name'))+
  geom_junction(data = ntrk2_junctions,
                aes(size = bdnf2_mean_psi),
                junction.y.max = 0.5)+
  geom_junction_label_repel(data = ntrk2_junctions,
                            aes(label = round(bdnf2_mean_psi, 2)),
                            junction.y.max = 0.5)+
  scale_size_continuous(range = c(0.1,1))+
  coord_cartesian(xlim = c(84668000, 84705000))

# 053--changes with bdnf (non-significant)

ntrk2_053_exons <- ntrk2_exons %>% 
  filter(transcript_id == 'ENST00000304053')
ntrk2_053_cds <- ntrk2_cds %>% 
  filter(transcript_id == 'ENST00000304053')

ntrk2_junctions <- splicing_junctions_2hr %>% 
  filter(gene_name == 'NTRK2') %>% 
  mutate(transcript_id = 'ENST00000304053')

ntrk2_053_exons %>% 
  ggplot(aes(xstart = start,
             xend= end,
             y = transcript_id))+
  geom_range(data = ntrk2_053_cds, aes(fill = type))+
  geom_range(fill = 'white',
             height = 0.25)+
  geom_range(data = ntrk2_053_cds)+
  geom_intron(data = to_intron(ntrk2_053_exons, 'transcript_name'))+
  geom_junction(data = ntrk2_junctions,
                aes(size = bdnf2_mean_psi),
                junction.y.max = 0.5)+
  geom_junction_label_repel(data = ntrk2_junctions,
                            aes(label = round(bdnf2_mean_psi, 2)),
                            junction.y.max = 0.5)+
  scale_size_continuous(range = c(0.1,1))+
  coord_cartesian(xlim = c(84668000, 84705000))
