#This is a tutorial to learn how to use the PCA function in R
#Source: StatQuest: PCA in R

#Creating a exemplary database

data.matrix <- matrix(nrow=100, ncol = 10)

colnames(data.matrix) <- c(
  paste("wt", 1:5, sep=""),
  paste("ko", 1:5, sep=""))

rownames(data.matrix) <- c(
  paste("gene", 1:100, sep=""))

#assigning random values

for (i in 1:100) {
  wt.values <- rpois(5, lambda = sample (x=10:1000, size = 1))
  ko.values <- rpois(5, lambda = sample (x=10:1000, size = 1))
  
  data.matrix [i,] <- c(wt.values, ko.values)
}

#Now we can perform PCA
#Note! The prcomp() function expects the sample to be rows and the genes to be columns, 
#therefore we need to tansponse the matrix

pca <- prcomp(t(data.matrix), scale = TRUE)

#pca contains three informations: x, sdev and rotation
#x contains the principal componetnts

plot(pca$x[,1], pca$x[,2])

#sdev contains how much standard deviation each component accounts for
#sdev^2 can be used to calculate how much variance each component accounts for
#Finally, we can transform that value into a percentage to better compare the PCs

pca.var <- pca$sdev^2
pca.var.per <- round(pca.var/sum(pca.var)*100, 1)

barplot(pca.var.per, main = "Scree Plot", xlab= "Principal Component", ylab= "Percentage of variation")

library(ggplot2)

pca.data <- data.frame(Sample=rownames(pca$x), X=pca$x[,1], Y=pca$x[,2] )

ggplot(data=pca.data, aes(x=X, y=Y, label=Sample)) +
 geom_text() +
 xlab(paste("PC1 -", pca.var.per[1], "%", sep =""))  +
 ylab(paste("PC2 -", pca.var.per[2], "%", sep =""))  +
 theme_bw()  + 
 ggtitle("My PCA Graph") 

#rotation shows how much influence(loading scores) each gene has on each PC. Positive influences push the points towards positive values and
#negative influences push the data points toward negative values. We'll now explore the genes' influence on PC1. Since the we're interested in the genes with the largest 
#influences we should sort them with the abs() function instead of sort them by value.

loading_scores <- pca$rotation[,1]
gene_scores <- abs(loading_scores)
gene_score_ranked <- sort(gene_scores, decreasing = TRUE)
top_10_genes <- names(gene_score_ranked[1:10])

top_10_genes

pca$rotation[top_10_genes,1]