#clustering
library(tidyverse)
library(cluster)
library(factoextra)
library(gridExtra)

new_ratio_bayesian_p_de <- read_csv("data/new_ratio_bayesian_p_de.csv")
new_ratio_bayesian_p_de <- new_ratio_bayesian_p_de %>% 
  mutate(gene = gsub("\\..*", "", gene)) 

never_sig <- new_ratio_bayesian_p_de %>% group_by(gene) %>% 
  filter(total_rna_sig == "not_significant") %>% count() %>% 
  filter(n == 3) %>% pull(gene)

log2clust = new_ratio_bayesian_p_de %>% 
  filter(!(gene %in% never_sig)) %>% 
  filter(!(is.na(gene_name))) %>% 
  # mutate(log2FoldNTR = log2(mean_bdnf_ntr / mean_control_ntr)) %>% 
  # filter(n_samp_passing_bdnf >2 & n_samp_passing_control > 2) %>% 
  dplyr::select(time,gene,gene_name,log2FoldChange) %>% 
  pivot_wider(names_from = 'time',
              values_from = c('log2FoldChange'), 
              values_fill = NA) %>% 
  drop_na() %>% 
  # mutate(
  #        across(where(is.numeric), ~ as.numeric(scale(.)))
  # ) %>%
  janitor::clean_names()

log2clust <- subset(log2clust, select = -c(gene))

log2clust <- log2clust %>% column_to_rownames('gene_name')
log2clust <- scale(log2clust)
#starting with the basics

kmeans2 <- kmeans(na.omit(log2clust), centers = 2, nstart = 25)
plot1 <- fviz_cluster(kmeans2, geom = 'point', data = log2clust)+
  ggtitle('k= 2')

kmeans3 <- kmeans(log2clust, centers = 3, nstart = 25)
plot2 <- fviz_cluster(kmeans3, geom = 'point',data = log2clust)+
  ggtitle('k = 3')

kmeans4 <- kmeans(log2clust, centers =4 , nstart = 25)
plot3<- fviz_cluster(kmeans4,geom = 'point', data = log2clust)+
  ggtitle('k = 4')

kmeans5 <- kmeans(log2clust, centers = 5, nstart = 25)
plot4 <- fviz_cluster(kmeans5, geom= 'point',data = log2clust)+
  ggtitle('k = 5')

grid.arrange(plot1,plot2,plot3,plot4, nrow = 2)

#elbow method for choosing clusters

fviz_nbclust(log2clust, kmeans, method = 'wss')

#shows i should use 7 clusters

#calculating gap statistic as another way to determine number of clusters
gap_stat <- clusGap(log2clust, FUN = kmeans, nstart = 25, K.max = 10, B = 50)
fviz_gap_stat(gap_stat)

#suggests 2 clusters should be used

final <- kmeans(log2clust, 7, nstart = 25)
fviz_cluster(final, geom = 'point', data = log2clust)

final2 <- kmeans(log2clust, 2, nstart = 25)
fviz_cluster(final2, data = log2clust)

#now im going to do the time series clustering

library(dtwclust)

part_clust <- tsclust(log2clust, type = 'partitional', k = 7L,
                      distance = 'dtw_basic', centroid = 'pam')
plot(part_clust)
