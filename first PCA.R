bdnf_full_counts_1_ <- read_excel("C:/Users/mlomb/OneDrive/Desktop/MRes project/bdnf_full_counts (1).xlsx")
bdnf_full_counts_1_ <- column_to_rownames(bdnf_full_counts_1_, "ensgene")
pca <- prcomp(t(bdnf_full_counts_1_),scale=TRUE)
plot(pca$x[,1], pca$x[,2])
pca.var <- pca$sdev^2
pca.var.per <- round(pca.var/sum(pca.var)*100,1)
barplot(pca.var.per, main="Scree Plot", xlab= "Principal Component", ylab= "Percent Variation")
library(ggplot2)
pca.data <- data.frame(Sample=rownames(pca$x), X=pca$x[,1], Y=pca$x[,2])
ggplot(data = pca.data, aes(x=X, y=Y, label=Sample)) +
  geom_text() +
  xlab(paste("PC1 - ", pca.var.per[1], "%", sep="")) +
  ylab(paste("PC2 - ", pca.var.per[2], "%", sep= "")) +
  theme_bw() +
  ggtitle("BDNF PCA Graph")
loading_scores <- pca$rotation[,1]
gene_scores <- abs(loading_scores)
gene_score_ranked <- sort(gene_scores, decreasing = TRUE)
top_100_genes <- names(gene_score_ranked[1:100])
pca$rotation[top_100_genes,1]