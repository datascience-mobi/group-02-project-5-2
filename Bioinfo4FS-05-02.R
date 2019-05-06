# Import dataset 
data <- readRDS(file = "C:/Users/alvar/Downloads/AML_gran_list.RDS")
annotation <- read.csv(file ="C:/Users/alvar/Downloads/sample_annotation.csv")

# Separate the dataframes into specific subgroups.
#   1st get the beta and the coverage values of the genes in 2 separate matrix
genes_beta <- data[["genes"]][,c(11:28)]
genes_coverage <- data[["genes"]][,c(29:46)]

# Let's try to do some QC with our dataset. For that, we'll make a distribution
# of the mean coverage of every gene. Genes with a very low coverage shall be removed. 
# Criteria for removal could be the lowest 10-percentile (?)
mean_coverage <- rowMeans(genes_coverage)
log_mean_coverage <- log10(mean_coverage)
hist(log_mean_coverage, breaks = 500, xlim = c(0,6))
abline(v = summary(log_mean_coverage)[2:5], col = c("blue", "red", "black", "orange"), lty = 2)
# This histogram shows the distribution of coverage values. 

# Let's try now to get the density function of methylation in one sample, 
# for example the first one (AML_myel_BM_S004XMA1_ctr.bed)
sample_01_beta <- genes_beta$AML_myel_BM_S004XMA1_ctr.bed 
sample_01_beta_density <- density(sample_01_beta, na.rm = TRUE) #Important: Removal of missing values 
plot(sample_01_beta_density)
# Histogram to compare with actual counts instead of relative density
hist(sample_01_beta, breaks = 100)

#Next we would plot all density functions of a group (e.g. cancer patients) together with the par command
# For cancer patients
par(mfrow=c(3,3))
samplename <- colnames(genes_beta)
for(i in 1:9){
  beta_dens <- density(genes_beta[,i], na.rm=T)
  plot(beta_dens, main = samplename[i], ylim=c(0,8))
}

# For the healthy patients
par(mfrow=c(3,3))
samplename <- colnames(genes_beta)
for(i in 10:18){
  beta_dens <- density(genes_beta[,i], na.rm=T)
  plot(beta_dens, main = samplename[i], ylim=c(0,8))
}

#Roughly, we can observe a higher density of methylation in cancer samples (mind the scale)
#How could we make a more quantitative comparison now? 

