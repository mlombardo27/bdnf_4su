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

#pulling data for GO

data_list <- part_clust@datalist
data_list <- as.data.frame(data_list)

col_names <- colnames(data_list)

part_clus_info <- tibble(col_names, clusters)

part_clus_info %>% 
  filter(clusters == '7') %>% 
  #filter(col_names == c('ARC','EGR1','FOS')) 
  pull(col_names) %>% 
  unique() %>% 
  clipr::write_clip()

#graphing GO data from string 

#part2

part_2_bp <- janitor::clean_names(read_tsv('part2_BP.tsv'))

part_2_bp %>% slice_max(strength, n = 10) %>% 
  mutate(term_description = fct_reorder(term_description,strength)) %>% 
  ggplot(aes(y = strength, fill = false_discovery_rate, x = term_description,size = observed_gene_count))+
  geom_point(pch = 21) +
  scale_fill_gradient(low = 'blue', high = 'red')+
  coord_flip()+ 
  ggpubr::theme_pubr()+
  theme(legend.position = 'right')+
  ggtitle("biological process GO on time series cluster k2")

part_2_mf <- janitor::clean_names(read_tsv('part2_MF.tsv'))

part_2_mf %>% slice_max(strength, n = 10) %>% 
  mutate(term_description = fct_reorder(term_description,strength)) %>% 
  ggplot(aes(y = strength, fill = false_discovery_rate, x = term_description,size = observed_gene_count))+
  geom_point(pch = 21) +
  scale_fill_gradient(low = 'blue', high = 'red')+
  coord_flip()+ 
  ggpubr::theme_pubr()+
  theme(legend.position = 'right')+
  ggtitle("molecular function GO on time series cluster k2")

part_2_cc <- janitor::clean_names(read_tsv('part2_CC.tsv'))

part_2_cc %>% slice_max(strength, n = 10) %>% 
  mutate(term_description = fct_reorder(term_description,strength)) %>% 
  ggplot(aes(y = strength, fill = false_discovery_rate, x = term_description,size = observed_gene_count))+
  geom_point(pch = 21) +
  scale_fill_gradient(low = 'blue', high = 'red')+
  coord_flip()+ 
  ggpubr::theme_pubr()+
  theme(legend.position = 'right')+
  ggtitle("cellular component GO on time series cluster k2")

#part5

part_5_bp <- janitor::clean_names(read_tsv('part5_BP.tsv'))

part_5_bp %>% slice_max(strength, n = 10) %>% 
  mutate(term_description = fct_reorder(term_description,strength)) %>% 
  ggplot(aes(y = strength, fill = false_discovery_rate, x = term_description,size = observed_gene_count))+
  geom_point(pch = 21) +
  scale_fill_gradient(low = 'blue', high = 'red')+
  coord_flip()+ 
  ggpubr::theme_pubr()+
  theme(legend.position = 'right')+
  ggtitle("biological process GO on time series cluster k5")

part_5_mf <- janitor::clean_names(read_tsv('part5_MF.tsv'))

part_5_mf %>% slice_max(strength, n = 10) %>% 
  mutate(term_description = fct_reorder(term_description,strength)) %>% 
  ggplot(aes(y = strength, fill = false_discovery_rate, x = term_description,size = observed_gene_count))+
  geom_point(pch = 21) +
  scale_fill_gradient(low = 'blue', high = 'red')+
  coord_flip()+ 
  ggpubr::theme_pubr()+
  theme(legend.position = 'right')+
  ggtitle("molecular function GO on time series cluster k5")

part_5_cc <- janitor::clean_names(read_tsv('part5_CC.tsv'))

part_5_cc %>% slice_max(strength, n = 10) %>% 
  mutate(term_description = fct_reorder(term_description,strength)) %>% 
  ggplot(aes(y = strength, fill = false_discovery_rate, x = term_description,size = observed_gene_count))+
  geom_point(pch = 21) +
  scale_fill_gradient(low = 'blue', high = 'red')+
  coord_flip()+ 
  ggpubr::theme_pubr()+
  theme(legend.position = 'right')+
  ggtitle("cellular component GO on time series cluster k5")


#part7

part_7_bp <- janitor::clean_names(read_tsv('part7_BP.tsv'))

part_7_bp %>% slice_max(strength, n = 10) %>% 
  mutate(term_description = fct_reorder(term_description,strength)) %>% 
  ggplot(aes(y = strength, fill = false_discovery_rate, x = term_description,size = observed_gene_count))+
  geom_point(pch = 21) +
  scale_fill_gradient(low = 'blue', high = 'red')+
  coord_flip()+ 
  ggpubr::theme_pubr()+
  theme(legend.position = 'right')+
  ggtitle("biological process GO on time series cluster k7")

part_7_mf <- janitor::clean_names(read_tsv('part7_MF.tsv'))

part_7_mf %>% slice_max(strength, n = 10) %>% 
  mutate(term_description = fct_reorder(term_description,strength)) %>% 
  ggplot(aes(y = strength, fill = false_discovery_rate, x = term_description,size = observed_gene_count))+
  geom_point(pch = 21) +
  scale_fill_gradient(low = 'blue', high = 'red')+
  coord_flip()+ 
  ggpubr::theme_pubr()+
  theme(legend.position = 'right')+
  ggtitle("molecular function GO on time series cluster k7")

part_7_cc <- janitor::clean_names(read_tsv('part7_CC.tsv'))

part_7_cc %>% slice_max(strength, n = 10) %>% 
  mutate(term_description = fct_reorder(term_description,strength)) %>% 
  ggplot(aes(y = strength, fill = false_discovery_rate, x = term_description,size = observed_gene_count))+
  geom_point(pch = 21) +
  scale_fill_gradient(low = 'blue', high = 'red')+
  coord_flip()+ 
  ggpubr::theme_pubr()+
  theme(legend.position = 'right')+
  ggtitle("cellular component GO on time series cluster k7")
#dbscan

library(dbscan)

# to plot the eps values
eps_plot <- kNNdistplot(log2clust, k=3)

# to draw an optimum line
eps_plot %>% abline(h = 0.6, lty = 2)

db <- dbscan(log2clust, eps = 0.8, minPts = 2)

fviz_cluster(db, log2clust, geom = 'point')

# k-medoid partitioning
# less sensitive to noise and outliers because it uses medoids as cluster centers instead of means 
#	“cluster”, “factoextra”

kmed_7 <- pam(log2clust, k = 7)
fviz_cluster(kmed_7, data = log2clust, geom= 'point')


kmed_3 <- pam(log2clust, k = 3)
fviz_cluster(kmed_3, data = log2clust, geom = 'point')

#getting data for GO

kmed_info <- as.data.frame(kmed_3$clustering) %>% 
  rownames_to_column('gene_name')

kmed_info %>% 
  filter(`kmed_3$clustering` == '3') %>% 
  #filter(gene_name == c('FOS','ARC','EGR1')) 
  pull(gene_name) %>% 
  unique() %>% 
  clipr::write_clip()

#graphing GO from string

#cluster3

kmed3_3_bp <- janitor::clean_names(read_tsv('kmed3_BP.tsv'))

kmed3_3_bp %>% slice_max(strength, n = 10) %>% 
  mutate(term_description = fct_reorder(term_description,strength)) %>% 
  ggplot(aes(y = strength, fill = false_discovery_rate, x = term_description,size = observed_gene_count))+
  geom_point(pch = 21) +
  scale_fill_gradient(low = 'blue', high = 'red')+
  coord_flip()+ 
  ggpubr::theme_pubr()+
  theme(legend.position = 'right')+
  ggtitle("biological process GO on time series cluster kmed cluster 3")

kmed3_3_mf <- janitor::clean_names(read_tsv('kmed3_MF.tsv'))

kmed3_3_mf %>% slice_max(strength, n = 10) %>% 
  mutate(term_description = fct_reorder(term_description,strength)) %>% 
  ggplot(aes(y = strength, fill = false_discovery_rate, x = term_description,size = observed_gene_count))+
  geom_point(pch = 21) +
  scale_fill_gradient(low = 'blue', high = 'red')+
  coord_flip()+ 
  ggpubr::theme_pubr()+
  theme(legend.position = 'right')+
  ggtitle("molecular function GO on time series cluster kmed cluster 3")

kmed3_3_cc <- janitor::clean_names(read_tsv('kmed3_CC.tsv'))

kmed3_3_cc %>% slice_max(strength, n = 10) %>% 
  mutate(term_description = fct_reorder(term_description,strength)) %>% 
  ggplot(aes(y = strength, fill = false_discovery_rate, x = term_description,size = observed_gene_count))+
  geom_point(pch = 21) +
  scale_fill_gradient(low = 'blue', high = 'red')+
  coord_flip()+ 
  ggpubr::theme_pubr()+
  theme(legend.position = 'right')+
  ggtitle("cellular component GO on time series cluster kmed cluster 3")


dyn_gene_kmed3<- kmed3_3_mf %>% 
  filter(grepl('dyn', term_description,ignore.case = TRUE))  %>% 
  pull(matching_proteins_in_your_network_labels) %>% 
  str_split(.,',')

dyn_gene_kmed3<- data.frame(gene_names = unlist(dyn_gene_kmed3))

new_ratio_bayesian_p_de %>% 
  filter(gene_name %in% dyn_gene_kmed3$gene_names) %>% 
  ggplot(aes(x = time, y= log2FoldChange, color= gene_name)) +
  geom_hline(yintercept = 0,size = 1)+
  geom_path(aes(x = time, y = log2FoldChange, group = gene_name), size = 3)+
  geom_point(size = 3,pch = 21)+
  scale_color_manual(values =randomcoloR::distinctColorPalette(k = 19))+
  ggtitle('log2foldChange for dynein related genes over time') +
  scale_x_continuous(labels = 0:6, breaks = 0:6)+
  labs(x = "Treatment Time")+
  ggpubr::theme_pubr()

dyn_gene_kmed3<- kmed3_3_cc %>% 
  filter(grepl('dyn', term_description,ignore.case = TRUE))  %>% 
  pull(matching_proteins_in_your_network_labels) %>% 
  str_split(.,',')

dyn_gene_kmed3<- data.frame(gene_names = unlist(dyn_gene_kmed3))

new_ratio_bayesian_p_de %>% 
  filter(gene_name %in% dyn_gene_kmed3$gene_names) %>% 
  ggplot(aes(x = time, y= log2FoldChange, color= gene_name)) +
  geom_hline(yintercept = 0,size = 1)+
  geom_path(aes(x = time, y = log2FoldChange, group = gene_name))+
  geom_point(size = 3,pch = 21)+
  ggtitle('log2foldChange for dynein related genes over time') +
  scale_x_continuous(labels = 0:6, breaks = 0:6)+
  labs(x = "Treatment Time")+
  ggpubr::theme_pubr()

#model based clustering

library(mclust)

#bic = bayesian information criterion

bic <- mclustBIC(log2clust)
plot(bic)

summary(bic)

#Gaussian finite mixture model fitted by EM algorithm 

mod1 <- Mclust(log2clust, x = bic)
plot(mod1, what = "classification")
summary(mod1, parameters = TRUE)

plot(mod1, what = 'uncertainty')

#icl = integrated completed likelihood, to select relevant numbers of classes
ICL <- mclustICL(log2clust)
summary(ICL)
plot(ICL)

#LRT = liklihood ratio test
#bootstrap sequential LRT for the number of mixture components
LRT <- mclustBootstrapLRT(log2clust, modelName = 'VVV')

#im going to initialize by using hierarchical clustering to get the partitions for the EM algorithm

#vvv model: ellipsoidal, varying volume, shape4 and orientation
#SVD: scaled SVD transformation
hc1 <- hc(log2clust, modelName = "VVV", use = "SVD")

bic1 <- mclustBIC(log2clust, initialization = list(hcPairs = hc1))

#vars: original variables
hc2 <- hc(log2clust, modelName = 'VVV', use = 'VARS')
bic2 <- mclustBIC(log2clust, initialization = list(hcPairs = hc2))

#eee model: ellipsoidal, equal volume, shape, and orientation
hc3 <- hc(log2clust, modelName = 'EEE', use = 'SVD')
bic3<- mclustBIC(log2clust, initialization = list(hcPairs = hc3))

#next we merge the BICs 

bic <- mclustBICupdate(bic1, bic2, bic3)
plot(bic)

mod1_update <- Mclust(log2clust, x = bic)
plot(mod1_update, what = 'classification')

# density estimation via gaussian finite mixture modelling

mod_density <- densityMclust(log2clust)
plot(mod_density, what = 'BIC')
plot(mod_denisty, what = 'density' , data = log2clust, breaks = 15)