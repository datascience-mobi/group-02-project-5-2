#Code 1: Bioinfo-05-02-QC_promoter.R

data <- AML_gran_list
annotation <- sample_annotation

#Split data in order to facilitate data cleanup
promoters <- data$promoters;
promoters_only <- promoters[,-c(1:10)]

beta <- promoters_only[,c(1:18)]
coverage<- promoters_only[,c(19:36)]

PatientInfo <- annotation[,c(3,4,11,29)]
rownames(PatientInfo) <- c(
  paste("AML", 1:9, sep=""),
  paste("CON", 1:9, sep="")
)



########################### DATA VISUALIZATION ######################

# In order to check for quality in the data, we need to set thresholds for minimum and 
# maximum coverage values. Too little coverage means unreliable beta values, whereas
# too high coverage could be caused by PCR duplicates

# An approach to setting the thresholds could be plotting the amount of NAs that we'd get
# for a set of threshold values and then choosing more or less arbitarily a proper threshold

# First let's visualize the general distribution of coverage values
par(mfrow=c(1,1))
coverage_tot_mean <- rowMeans(coverage) 
log_coverage_mean <- log10(coverage_tot_mean)

hist(log_coverage_mean, breaks = 250, xlim = c(0,5.5))
abline(v=quantile(log_coverage_mean, probs = c(seq(0.1, 1.0, by= 0.1))), col="red")
abline(v=quantile(log_coverage_mean, probs = c(0.05, 0.95), col="cyan"), lty=2)
# The green lines tells us that the 5% quantile is around cov=100 
# and the 95% quantile aorund cov=5.500

quantile(coverage_tot_mean, probs=.05)
#Result: 99.33333
quantile(coverage_tot_mean, probs=.95)
#Result: 5620 


########################### QUALITY CONTROL #############################

# Chromosomes X & Y have to be removed from the analysis
cover_nox <- coverage[-which(promoters$Chromosome == "chrX"), ]
cover_noxy <- cover_nox[-which(promoters$Chromosome == "chrY"), ]
rm(coverage, cover_nox)
beta_nox <- beta[-which(promoters$Chromosome == "chrX"), ]
beta_noxy <- beta_nox[-which(promoters$Chromosome == "chrY"), ]
rm(beta, beta_nox)

print( (nrow(beta_noxy)*100) / nrow(promoters) ) ##RESULT = 95.93893%


################################  UPPER threshold  ########################################

# In order to determine where to set a threshold, we create a plot that tells us
# how many values we lose depending on where we set the threshold

coverage_histo <- unlist(cover_noxy) #This allows to inspect the whole dataframe freely

NAs_up <- c()
thresholds_up <- c()
for (threshold in seq(0, max(coverage_histo), by = 200)) {
  thresholds_up <- append(thresholds_up, threshold)
  Unreliable <- coverage_histo[which(coverage_histo > threshold)]
  N_Unreliable <- length(Unreliable)
  NAs_up <- append(NAs_up, N_Unreliable)
}

length(NAs_up)
length(thresholds_up) #Checking if both have the same length so that they can be plotted

# Let's try making the x axis log10 to stretch the right part

par(mfrow = c(1, 1))
plot(
  x = log10(thresholds_up),
  y = NAs_up,
  type = "l",
  main = "Amount of Beta values converted to NA by threshold",
  xlab = "threshold position",
  ylab = "Unreliable Values"
)
abline(v = quantile(log_coverage_mean, probs = c(seq(0.95, 1.0, by = 0.01))), 
       col = "green"
)
abline(
  v = quantile(log_coverage_mean, probs = c(seq(0, 1.0, by = 0.1))),
  col = "red",  lty = 3
)

# PROBLEM: we are not taking into account the amount of NAs that were already in the Beta values as NaN

beta_NaN <- promoters_only[which(promoters_only[,c(1:18)]== "NaN"),] #Select the genes that have at least 1 NaN

colorder <- c() 
for(i in 1:18){
  colorder <- append(colorder, i)
  k <- i+18
  colorder <- append(colorder, k)
}
#This little for loop generates the vector of order that will sort the data frame putting the 
#beta value next to the coverage value of each sample instead of all betas and then all coverages
beta_NaN_table <- beta_NaN[, colorder] #Now the data frame has the beta values next to the coverage
#We observe that NaN values correspond to a coverage of 0.

is.na(beta_NaN_table[c(1:20), 1])
sum(is.na(beta_NaN_table[c(1:20), 1])) #NaN are treated as NAs by the is.na function

NAs_given <- c()
for(i in 1:18){
  NAs_given <- append( NAs_given, sum( is.na (promoters_only[, i]) ) )  
  
}
sum(NAs_given)

percent_of_NAs_given <- sum(NAs_given)/(nrow(promoters_only)*18)

##################################################################################

# On a second thought, it would be more meaningful to check how many 
# ROWS (=promotors) we lose depending on where we set the threshold

############### LOWER THRESHOLD 

# (Ziller, M. J., Hansen, K. D., Meissner, A., & Aryee, M. J. (2015)).
# Coverage recommendations for methylation analysis by whole-genome bisulfite sequencing.
# Nature methods, 12(3), 230.)


lower_threshold <- 30

beta_low_trimmed <- beta_noxy
for (j in 1:ncol(beta_low_trimmed)) {
  for (i in 1:nrow(beta_low_trimmed)) {
    if (cover_noxy[i, j] < lower_threshold) {
      beta_low_trimmed[i, j] <- NA
    }
  }
}

beta_low_NA_rm <- beta_low_trimmed[-which(rowSums(is.na(beta_low_trimmed)) > 2), ]
which(rowSums(is.na(beta_low_NA_rm)) > 2) ## RESULT = 0

print( (nrow(beta_low_NA_rm)*100) / nrow(beta_noxy))
print( (nrow(beta_low_NA_rm)*100) / nrow(promoters_only)) ##RESULT=91.02

############### UPPER THRESHOLD 

coverage_histo <- unlist(cover_noxy) #This allows to inspect the whole dataframe freely

rows_lost <- c()
rows_lost_tot <- c(0)
thresholds_2NAs <- c()
df <- cover_noxy
for(threshold in quantile(coverage_histo, probs=seq(1, 0.9, by=-0.005))){
  rows_lost <- c()
  for(i in 1:nrow(df) ){
    
    AML_position_above_limit <- which(df[i,c(1:9)] > threshold) #Gives the position of the AML values above threshold
    CON_position_above_limit <- which(df[i,c(10:18)] > threshold) #Same as above but with control samples
    
    if( length(AML_position_above_limit > 2) | length(CON_position_above_limit > 2) ) {
      
      rows_lost <- append(rows_lost, i)
    }
  }
  
  rows_lost_tot <- append( rows_lost_tot, ( rows_lost_tot[length(rows_lost_tot)] + length(rows_lost) ) )
  thresholds_2NAs <- append(thresholds_2NAs, threshold)
  print(rows_lost_tot)
  print(thresholds_2NAs)
  if(length(rows_lost) > 0){
    df <- df[-rows_lost,]
  }
  print(nrow(df))
}
rows_lost_tot <- rows_lost_tot[-1]

#Rows lost vs threshold
plot(x=thresholds_2NAs, y=rows_lost_tot, type="b") 

#Percentage of information lost
rows_lost_per <- (rows_lost_tot*100) / nrow(cover_noxy)
plot(x=thresholds_2NAs, y=rows_lost_per, type="b") 


#Rows lost vs threshold, better visualization without last value
plot(x=thresholds_2NAs[2:21], y=rows_lost_tot[2:21], type="b")

#Percentage of information lost, better visualization without last value
rows_lost_per <- (rows_lost_tot*100) / nrow(cover_noxy)
plot(x=thresholds_2NAs[c(2:21)], y=rows_lost_per[c(2:21)], type="b") 
abline(h=6, lty=2)

#UP until the removal of coverage values below 30, we've lost 10% of the information.
#Accordingly, the upper threshold should be chosen so that no more than 15% of the information (rows)
#get lost.

##  If we follow that, the upper threshold would be the 99.5% quantile, in which we would lose 4.56% of the information,
# leaving us with 86.5%

################################  APPLYING THRESHOLDS TO BETA VALUES  ########################################
# Now that we've decided where to set our coverage thresholds, every beta value whose
# coverage is either too low or too high shall be converted to NA

lower_threshold <- 30
upper_threshold <- quantile(coverage_histo, probs=0.995)

beta_trimmed <- beta_noxy
for (j in 1:ncol(beta_trimmed)) {
  for (i in 1:nrow(beta_trimmed)) {
    if (cover_noxy[i, j] < lower_threshold |
        cover_noxy[i, j] > upper_threshold) {
      beta_trimmed[i, j] <- NA
    }
  }
}
View(beta_trimmed)


# Removing genes that have more than 2 NAs in individual cohorts
beta_AML <- beta_trimmed[,c(1:9)]
beta_con <- beta_trimmed[,c(10:18)]

NAs_AML <- rowSums(is.na(beta_AML))
NAs_con <- rowSums(is.na(beta_con))

AML <- cbind(beta_AML, NAs_AML)
con <- cbind(beta_con, NAs_con)

beta_qc <- cbind(AML, con)
beta_partly_cleaned <- beta_qc[-which(beta_qc[, 10] > 2 | beta_qc[, 20] > 2), ]

# Replacing NAs with the mean value of the cohort

for (i in 1:nrow(beta_partly_cleaned)) {
  n <- beta_partly_cleaned[i,10]
  
  if (n==1 | n==2) {
    for (j in 1:9){
      if (is.na(beta_partly_cleaned[i,j])){
        beta_partly_cleaned[i,j] <- rowMeans(beta_partly_cleaned[i,1:9], na.rm = T)
      }
    }
  }
  
  n <- beta_partly_cleaned[i,20]
  
  if (n==1 | n==2) {
    for (j in 11:19){
      if (is.na(beta_partly_cleaned[i,j])){
        beta_partly_cleaned[i,j] <- rowMeans(beta_partly_cleaned[i,11:19], na.rm = T)
      }
    }
  }
}

beta_both_clean <- beta_partly_cleaned[,-c(10,20)]


### Data Normalization 

for (j in 1:18) {
  for (i in 1:nrow(beta_both_clean)) {
    beta_both_clean[[i,j]] <- log2(beta_both_clean[i,j] / (1-beta_both_clean[i,j]) )
    
  }
}

promotors_normalized <- beta_both_clean

View(promotors_normalized)

write.table(promotors_normalized, file="promoters_normalized.csv",sep=",",row.names=T) #Export table
promoters_norm <- read.csv(file="promoters_normalized.csv")




promotors_normalized <- read.csv(file = "promoters_normalized.csv")

######################PCA#############################

#The data must be turned into a matrx to use the function prcomp (Principal components)

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

############### Tests #################

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




############ Feature selction PC2 ####################

loading_scores <- promotors_pca$rotation
gene_scores <- abs(loading_scores)

PC2gene_scoresranked <- sort(gene_scores[,2], decreasing = TRUE)

par(mfrow=c(1,3)) #Overview of PC2 
plot(PC2gene_scoresranked)
abline(v=c(0,10000), lty=2)
plot(PC2gene_scoresranked[1:10000])
abline(v=c(0), col="green", lty=2)
plot(PC2gene_scoresranked[1:3000])

dev.off()

plot(PC2gene_scoresranked[1:5000])
abline(v=c(200, 500), col="green", lty=2)



GOI_200 <- PC2gene_scoresranked[1:200]
GOI_200_names <- names(GOI_200)

GOI_200_cluster <- promotors_matrix_noInf[GOI_200_names,]
GOI_200_k <- kmeans(t(GOI_200_cluster), centers=2)

View(GOI_200_k$cluster)


GOI_500 <- PC2gene_scoresranked[1:500]
GOI_500_names <- names(GOI_500)

GOI_500_cluster <- promotors_matrix_noInf[GOI_500_names,]
GOI_500_k <- kmeans(t(GOI_500_cluster), centers=2)

View(GOI_500_k$cluster)

GOI_Everything <- kmeans(t(promotors_matrix_noInf), centers=2)

View(GOI_Everything$cluster)


GOI_22750 <- PC2gene_scoresranked[1:22750]
GOI_22750_names <- names(GOI_22750)
GOI_22750_cluster <-promotors_matrix_noInf[GOI_22750_names,]

GOI_22750_k <- kmeans(t(GOI_22750_cluster), centers=2)

View(GOI_22750_k$cluster)

#Patient AML 3 is an outsider


GOI_FeatureSelection <- GOI_22750_cluster[,-3]

View(GOI_FeatureSelection)

GOI_FeatureSelection_k <- kmeans(t(GOI_FeatureSelection), centers=2)

View(GOI_FeatureSelection_k$cluster)

#Tadaaaaaaaaaaaaaaaaaah

#Now let's look for DMRs
#First we run a t test on the selected genes

Ttest_results <- data.frame(Stat = rep(0,22750),pval = rep(0,22750) , mean_difference = rep(0,22750))

for (i in 1:nrow(GOI_FeatureSelection)) {
  
  x = GOI_FeatureSelection[i,1:8]
  y= GOI_FeatureSelection[i,9:17]
  tt <- t.test(x,y)
  Ttest_results[i,1] <- tt$statistic
  Ttest_results[i,2] <- tt$p.value
  Ttest_results[i,3] <- mean(x) - mean(y)
  
}
rownames(Ttest_results) <- rownames(GOI_FeatureSelection)
 
View(Ttest_results)


plot(Ttest_results$mean_difference, -log10(Ttest_results$pval), main = "Volcano Plot: T-test Results")

#Let's put a threshold to the amount of difference in methylation that we will consider

plot(rowMeans(GOI_FeatureSelection[,9:17]), rowMeans(GOI_FeatureSelection[,1:8]))
abline(a=0, b=1)


Ttest_results_rankedMean <- Ttest_results[order(abs(Ttest_results$mean_difference)),]

plot(Ttest_results_rankedMean$mean_difference )

abline(h=c(1,1.5, 1.75))

plot(Ttest_results_rankedMean[20000:23000,])

abline(h=c(1, 1.5, 2, 2.5, 3, 3.5))

hist(Ttest_results_rankedMean$mean_difference)

Ttest_results_trimmedMean <- Ttest_results_rankedMean[which(Ttest_results_rankedMean$mean_difference >= 1|Ttest_results_rankedMean$mean_difference <= -1 ),]

length(t(Ttest_results_trimmedMean)) 



#Let's now examine the data by pvalue 

Ttest_results_rankedPValue <- Ttest_results_trimmedMean[order(Ttest_results_trimmedMean$pval),]

View(Ttest_results_rankedPValue)


hist(Ttest_results_rankedPValue$pval, main = "DMRs by p.value", xlab = "p.value")

hist(Ttest_results_rankedPValue$pval[1:10000], main = "DMRs by p.value", xlab = "p.value")


Ttest_results_rankedPValue_001 <- Ttest_results_rankedPValue[which(Ttest_results_rankedPValue$pval <= 0.01 ),]

Ttest_results_rankedPValue_005 <- Ttest_results_rankedPValue[which(Ttest_results_rankedPValue$pval <= 0.05 ),]

Ttest_results_rankedPValue_01 <-  Ttest_results_rankedPValue[which(Ttest_results_rankedPValue$pval <= 0.1  ),]

plot(x=c(0.01, 0.05, 0.1), c(length(t(Ttest_results_rankedPValue_001)), 
       length(t(Ttest_results_rankedPValue_005)), 
       length(t(Ttest_results_rankedPValue_01)) 
       ), 
     xlab = "P.Value", ylab = "Promotors_Left")



plot(Ttest_results$mean_difference, -log10(Ttest_results$pval), main = "Volcano Plot: T-test Results",  xlab = "Mean difference", ylab = "-log10(p.value)")

abline(v = 1.5, lty = 2, lwd = 3)
abline(v = -1.5, lty = 2, lwd = 3)
abline(h = -log10(0.05), lty = 2, lwd = 3)

points( Ttest_results_rankedPValue_01[,3], -log10(Ttest_results_rankedPValue_01[,2]), col = "red")
points( Ttest_results_rankedPValue_005[,3], -log10(Ttest_results_rankedPValue_005[,2]), col = "green")
points( Ttest_results_rankedPValue_001[,3], -log10(Ttest_results_rankedPValue_001[,2]), col = "blue")


###### Multiple comparison correction #######

#We're going to use the Holm-Šidák correction which is more powerful than the Holm-Bonferroni
# correction and the Bonferroni correction

#Since we have to repeat an algorhithm until we found out first acceptable null-hypothesis 
#we're going to use a while-loop

Fail_reject <- 0
Index <- 1
alpha <- 0.05
m <- length(t(Ttest_results_rankedPValue))
            
while (Fail_reject == 0) {
  
 
  alpha2 = 1-((1-alpha)^(1/(m-(Index-1))))
  
  if(Ttest_results_rankedPValue[Index,2] > alpha2 ){
    
    Fail_reject <- Index
    
    
  }else{
    
   Index = Index+1
  }
  
}

View(Fail_reject) #At 133 the p value start become unreliable

Reliable_Ttest_results_p <- Ttest_results_rankedPValue[1:132,]
Reliable_Ttest_results_fold <- Reliable_Ttest_results_p[order(Reliable_Ttest_results_p$mean_difference),]



#Manhattan plot


install.packages("qqman")
library(qqman)

# A Manhattan plot has basically 3 elements: Chromosome, Basepair and p value
# Because the BP and CHR are in the original promoters dataset before the quality control,
# we'll have to select only the rows that made it through the QC

promoters <- data$promoters;


names <- rownames(Reliable_Ttest_results) # This gives us a vector with the genes that survived the QC
promotersQC <- promoters[names, ] # We pick only the rows with those genes 

chr <- promotersQC$Chromosome
bp <- promotersQC$Start # The choice of Start was arbitrary. End could also be used, even the median of Start and End if you wanted. I don't think it changes much
p <-  Reliable_Ttest_results$pval

df <- as.data.frame(cbind(chr, bp, p)) # If we keep it as a dataframe, then later we can refer to the columns by writing their names 
rownames(df) <- names

manhattan(df, chr="chr", bp="bp", p="p", suggestiveline = F, genomewideline = -log10(Reliable_Ttest_results[20,2])) 

#We'll pick the top 20 genes

Selected_promoters_p <- promotersQC[rownames(Reliable_Ttest_results_p[1:20,]),]
Selected_promoters_fold <- promotersQC[rownames(Reliable_Ttest_results_fold[1:20,]),]

View(Selected_promoters_p)
View(Selected_promoters_fold)





