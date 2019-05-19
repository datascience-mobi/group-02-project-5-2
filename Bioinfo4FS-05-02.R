# Import dataset 
data <- readRDS(file = "C:/Users/alvar/Downloads/AML_gran_list.RDS")
annotation <- read.csv(file ="C:/Users/alvar/Downloads/sample_annotation.csv")

# Separate the dataframes into specific subgroups.
#   1st get the beta and the coverage values of the genes in 2 separate matrix
genes_beta <- data[["genes"]][,c(11:28)]
genes_coverage <- data[["genes"]][,c(29:46)]

##################################    QC    #########################################

##################################  [ COVERAGE

# Let's try to do some QC with our dataset. For that, we'll make a distribution
# of the mean coverage of every gene. Genes with a very low coverage shall be removed. 
# Criteria for removal shall be determined according to the distribution
par(mfrow=c(1,1))

mean_coverage <- rowMeans(genes_coverage) 
log_mean_coverage <- log10(mean_coverage)
sd_coverage <- apply(genes_coverage, 1, sd)
log_sd_coverage <- log10(sd_coverage)
sd_coverage_dens <- density(log_sd_coverage)

# Distribution looks bimodal, maybe it shouldn't?
# Trying to split the data into case and control and then visualize each distribution
coverage_mean_cancer <- rowMeans(genes_coverage[,c(1:9)])
coverage_mean_control <- rowMeans(genes_coverage[,c(10:18)]) 
par(mfrow=c(1,2)) 
hist(coverage_mean_cancer, breaks = 500, xlim = c(0,200000))
hist(coverage_mean_control, breaks = 500, xlim = c(0,200000))

log_cover_mean_cancer <- log10(coverage_mean_cancer)
log_cover_mean_control <- log10(coverage_mean_control)
par(mfrow=c(1,2)) 
hist(log_cover_mean_cancer, breaks = 200, xlim = c(0,6))
abline(v=quantile(log_cover_mean_cancer, probs = c(seq(0.1, 1.0, by= 0.1))), col="red")
hist(log_cover_mean_control, breaks = 200, xlim = c(0,6))
abline(v=quantile(log_cover_mean_control, probs = c(seq(0.1, 1.0, by= 0.1))), col="green")

#Still looks bimodal. Trying each patient in a cohort
#For cancer patients
par(mfrow=c(3,3))
samplename <- colnames(genes_coverage)
for(i in 1:9){
  log_cover_ca <- log10(genes_coverage[,i])
  hist(log_cover_ca, main = samplename[i], breaks=200)  #ylim=c(0,8)
}

# For control patients
par(mfrow=c(3,3))
samplename <- colnames(genes_coverage)
for(i in 10:18){
  log_cover_co <- log10(genes_coverage[,i])
  hist(log_cover_co, main = samplename[i], breaks=200)  #ylim=c(0,8)
}

#General comparison - too large :/
par(mfrow=c(6,3))
samplename <- colnames(genes_coverage)
for(i in 1:18){
  log_cover_all <- log10(genes_coverage[,i])
  hist(log_cover_all, main = samplename[i], breaks=200)  #ylim=c(0,8)
}
#Still looks kindof bimodal. Maybe it's not so wrong?




# I'm trying to include the standard deviation and the mean coverage values of each gene in a single plot, but I seem to be missing something.
hist(log_mean_coverage, breaks = 500, xlim = c(0,6))
abline(v = summary(log_mean_coverage)[2:5], col = c("blue", "red", "black", "orange"), lty = 2)
lines(log_sd_coverage, col="green")

hist(log_sd_coverage, breaks=200, xlim = c(0,6))
lines(sd_coverage_dens, col="red")

################################## COVERAGE ]


##################################  [ METHYLATION : BETA

# Let's try now to get the density function of methylation in one sample, 
# for example the first one (AML_myel_BM_S004XMA1_ctr.bed)
sample_01_beta <- genes_beta$AML_myel_BM_S004XMA1_ctr.bed 
sample_01_beta_density <- density(sample_01_beta, na.rm = TRUE) #Important: Removal of missing values 
plot(sample_01_beta_density)
# Histogram to compare with actual counts instead of relative density
hist(sample_01_beta, breaks = 100)
#Now both at the same time huehue
par(mfrow=c(1,2))
plot(sample_01_beta_density)
hist(sample_01_beta, breaks = 100)

#Next we would plot all density functions of a group together using the par command
# For cancer patients
par(mfrow=c(3,3))
samplename <- colnames(genes_beta)
for(i in 1:9){
  beta_dens <- density(genes_beta[,i], na.rm=T)
  plot(beta_dens, main = samplename[i], ylim=c(0,8))
}

# For control patients
par(mfrow=c(3,3))
samplename <- colnames(genes_beta)
for(i in 10:18){
  beta_dens <- density(genes_beta[,i], na.rm=T)
  plot(beta_dens, main = samplename[i], ylim=c(0,8))
}

#Roughly, we can observe a higher density of methylation in cancer samples 
#How could we make a more quantitative comparison now? 

################################## METHYLATION : BETA ]


##################################    NORMALIZATION    #########################################

# Conversion factor between beta and m values: M = log2(Beta/1-Beta)

# Let's use the sample patient (cancer) we selected before to see how normalization looks like
sample_01_m <- sample_01_beta
for (i in 1:56234) { 
  sample_01_m[[i]] <- log2(sample_01_beta[[i]]/(1-sample_01_beta[[i]]))
                      }
hist(sample_01_m, breaks = 100) 

# Let's try with a non cancer patient now (first patient of the control group)
sample_02_beta <- genes_beta$Neut_band_BM_S00JGXA1.bed
sample_02_m <- sample_02_beta
for (i in 1:56234) { 
  sample_02_m[[i]] <- log2(sample_02_beta[[i]]/(1-sample_02_beta[[i]]))
}
hist(sample_02_m, breaks = 100, xlim=c(-10,10)) 

# Now let's try to get the M values of the first gene accross all samples, 
# and then plot the distribution of each cohort (cancer and then control)
gene_01_beta <- genes_beta[1,]
gene_01_m <- gene_01_beta
for (i in 1:18) { 
  gene_01_m[[i]] <- log2(gene_01_beta[[i]]/(1-gene_01_beta[[i]])) # Transform all Beta into M values
}
gene_01_m_t <- as.data.frame(t(gene_01_m)) # Convert the row into a column so that histogram can be plotted
hist(gene_01_m_t$ENSG00000223972, breaks = 18) # Preliminary histogram of the distribution of the M value of all samples (case & control)
shapiro.test(gene_01_m_t$ENSG00000223972) # Run a Shpairo test to check for normality of the data

#Convert all Beta values into M values
genes_m <- genes_beta
for (j in 1:18) {
  for (i in 1:56234) {
    genes_m[[i,j]] <- log2(genes_beta[i,j] / (1-genes_beta[i,j]) )
    
  }
}
genes_m_can <- genes_m[,c(1:9)]
genes_m_can_mean <- rowMeans(genes_m_can)
genes_m_con <- genes_m[,c(10:18)]
genes_m_con_mean <- rowMeans(genes_m_con)
par(mfrow=c(1,2))
hist(genes_m_can_mean, breaks = 250, ylim=c(0,1500), xlim=c(-8,8), main="Cancer patients: Mean M values")
axis(side=1, at=seq(-8, 8, by=1), tick=T)
abline(v=2.5, col="red")
hist(genes_m_con_mean, breaks = 250, ylim=c(0,1500), xlim=c(-8,8), main="Healthy patients: Mean M values")
axis(side=1, at=seq(-5, 5, by=1), tick=T)
abline(v=3, col="orange")
