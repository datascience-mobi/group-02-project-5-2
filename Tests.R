promotors_normalized <- read.csv(file = "promoters_normalized.csv")

#PCA

#The data must be turned into a matrx to use the fucntion prcomp (Principal components)

promotors_matrix <- as.matrix(promotors_normalized) 

class(promotors_matrix)

#Beta values of 0 and 1 have been respectively turned into -Inf and Inf values, the PCA is not able to work with these values
#Ina and -Inf values are turned turned into extremely low/high values for M values, respectively -20 and 20

promotors_matrix_noInf <- promotors_matrix

promotors_matrix_noInf[which(promotors_matrix_noInf == Inf)] <- 20
promotors_matrix_noInf[which(promotors_matrix_noInf == -Inf)] <- -20

#In order to simplify our dataset the names of the sample are turned into intuitive names

colnames(promotors_matrix_noInf) <- c(
  paste("AML", 1:9, sep=""),
  paste("CON", 1:9, sep=""))

View(promotors_matrix_noInf)

#The PCA is run

promotors_pca <- prcomp(t(promotors_matrix_noInf), scale = TRUE)

sample_annotation <- read.csv(file="sample_annotation.csv")

#we create a dataframe with possible sources of variation and the PCs to inspect them closely 
# for now we want to inspect the effects that can be quantified with the kruskal-wallis test 

pca_promotors <- promotors_pca$x

PatientInfo_kruskal <- sample_annotation[,c(3,28,29,33)]

kruskal_df <-cbind(PatientInfo_kruskal,pca_promotors)
colnames(kruskal_df) <- c(1:22)


# we want to apply the kruskal wallis test on PC1-5(col 5-10 in the data frame)
# with the effects cellType, Disease, Donor-ID and Biomaterial_Provider (col 1-4 in data frame)


# we create a matrix with p.values with PCs as coloumns and the effects as rows
heatmap_kruskal <- matrix(nrow=4,ncol=5)

for (j in 1:4){
  for (i in 5:19) {
    heatmap_kruskal[j,i-4] <- kruskal.test(kruskal_df[,i] ~ kruskal_df[,j], data = kruskal_df)$p.value
    
  }
}
# now we create a similiar matrix with the statistic value

kruskal_statistic <- matrix(nrow=4,ncol=5)

for (j in 1:4){
  for (i in 5:19) {
    kruskal_statistic[j,i-4] <- kruskal.test(kruskal_df[,i] ~ kruskal_df[,j], data = kruskal_df)[["statistic"]][["Kruskal-Wallis chi-squared"]]
    
  }
}

# now we want to aplly the wilcoxon-rank sum test on the same PCs with the effect gender

gender <- sample_annotation[,36]
wilcoxon_df <- cbind(gender,pca_promotors)
colnames(Wilcoxon_df) <- c(1:19)

# with this code wa have a matrix with the p.values 
heatmap_wilcoxon <- matrix(nrow=1,ncol=5)

  for (i in 2:6) {
    heatmap_wilcoxon[1,i-1] <- wilcox.test(wilcoxon_df[,i] ~ wilcoxon_df[,1], data = wilcoxon_df)$p.value
    
  }
# and here the statistic_values

wilcoxon_statistic <- matrix(nrow=1,ncol=5)


  for (i in 2:6) {
    wilcoxon_statistic[1,i-1] <- wilcox.test(wilcoxon_df[,i] ~ wilcoxon_df[,1], data = wilcoxon_df)[["statistic"]][["W"]]
    
  }
  
  ############################### CORRELATION + PERMUTATION #############################

age <- sample_annotation[,c(34)]
age_separated <- strcapture("(.*)-(.*)", as.character(age), data.frame(type_1 = "", type_2 = ""))
age <- age_separated[,1]
correlation_df <- cbind(age, pca_promotors)
heatmap_correlation <- matrix(nrow=1,ncol=5)

#
# PERMUTATION TEST!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
#


cor.perm <- function (x,y, nperm = 499)
{
  r.obs <- cor (x = x, y = y)
  P.par <- cor.test (x = x, y = y)$p.value
  #  r.per <- replicate (nperm, expr = cor (x = x, y = sample (y)))
  r.per <- sapply (1:nperm, FUN = function (i) cor (x = x, y = sample (y)))
  r.per <- c(r.per, r.obs)
  hist (r.per, xlim = c(-1,1))
  abline (v = r.obs, col = 'red')
  P.per <- sum (abs (r.per) >= abs (r.obs))/(nperm + 1) 
  return (P.per = P.per)
}

for (i in 2:6){
  heatmap_correlation[1, i-1] <- cor.perm(x=correlation_df[,1], y=correlation_df[,i])
  

}
  
correlation_statistics <- matrix(nrow=1, ncol=5)

cor.perm_statistics <- function (x,y, nperm = 499)
{
  r.obs <- cor (x = x, y = y)
  P.par <- cor.test (x = x, y = y)[["statistic"]][["t"]]
  #  r.per <- replicate (nperm, expr = cor (x = x, y = sample (y)))
  r.per <- sapply (1:nperm, FUN = function (i) cor (x = x, y = sample (y)))
  r.per <- c(r.per, r.obs)
  hist (r.per, xlim = c(-1,1))
  abline (v = r.obs, col = 'red')
  P.per <- sum (abs (r.per) >= abs (r.obs))/(nperm + 1) 
  return (P.per = P.per)
}

for (i in 2:6){
  correlation_statistics[1, i-1] <- cor.perm_statistics(x=correlation_df[,1], y=correlation_df[,i])
  
  
}




#################################### HEATMAP ###################################


heatmap <- rbind(heatmap_kruskal,heatmap_wilcoxon,heatmap_correlation)
colnames(heatmap) <- c("PC1","PC2","PC3","PC4","PC5")
rownames(heatmap) <- c("cellType","Disease","Biomaterial-provider","Donor-ID","gender","Age")

statistics <- rbind(kruskal_statistic,wilcoxon_statistic,correlation_statistics)
colnames(statistics) <- c("PC1","PC2","PC3","PC4","PC5")
rownames(statistics) <- c("cellType","Disease","Biomaterial-provider","Donor-ID","gender","Age")


library(gplots)

my_palette <- colorRampPalette(c("red", "yellow"))(n = 3)
col_breaks <- c(seq(0,0.01,length=2),  
                seq(0.0105, 1,length=2)
) 

png("./heatmap_pvalues9.png",    # create PNG for the heat map        
    width = 5*300,        # 5 x 300 pixels
    height = 5*300,
    res = 300,            # 300 pixels per inch
    pointsize = 8)

par(mar=c(1,1,1,1))
heatmap.2(heatmap, 
          main = "Significance of effects on data",
          density.info="none",
          trace="none",
          lwid= c(2,5),
          margins =c(6,10),
          col=my_palette,
          breaks=col_breaks,
          dendrogram = "none" ,
          Colv = "NA"
)

dev.off()

my_palette <- colorRampPalette(c("red", "yellow"))(n = 3)
col_breaks <- c(seq(0,0.05,length=2),  
                seq(0.0505, 1,length=2)
) 

png("./heatmap_pvalues10.png",    # create PNG for the heat map        
    width = 5*300,        # 5 x 300 pixels
    height = 5*300,
    res = 300,            # 300 pixels per inch
    pointsize = 8)

par(mar=c(1,1,1,1))
heatmap.2(heatmap, 
          main = "Significance of effects on data",
          density.info="none",
          trace="none",
          lwid= c(2,5),
          margins =c(6,10),
          col=my_palette,
          breaks=col_breaks,
          dendrogram = "none" ,
          Colv = "NA"
)

dev.off()

my_palette <- colorRampPalette(c("red", "yellow"))(n = 3)
col_breaks <- c(seq(0,0.05,length=2),  
                seq(0.0505, 1,length=2)
) 

png("./heatmap_statistics2.png",    # create PNG for the heat map        
    width = 5*300,        # 5 x 300 pixels
    height = 5*300,
    res = 300,            # 300 pixels per inch
    pointsize = 8)

par(mar=c(1,1,1,1))
heatmap.2(statistics, 
          main = "Significance of effects on data",
          density.info="none",
          trace="none",
          lwid= c(2,5),
          margins =c(6,10),
          col=my_palette,
          breaks=col_breaks,
          dendrogram = "none" ,
          Colv = "NA"
)

dev.off()

