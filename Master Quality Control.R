#Code 1: Bioinfo-05-02-QC_promoter.R

data <- readRDS(file="AML_gran_list.RDS")
annotation <- read.csv(file="sample_annotation.csv")

#Split data in order to facilitate data cleanup
promoters <- data$promoters;
promoters_only <- promoters[,-c(1:10)]

beta <- promoters_only[,c(1:18)]
coverage<- promoters_only[,c(19:36)]

########################### QUALITY CONTROL - COVERAGE ######################

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

# Chromosomes X & Y have to be removed from the analysis
cover_nox <- coverage[-which(promoters$Chromosome == "chrX"), ]
cover_noxy <- cover_nox[-which(promoters$Chromosome == "chrY"), ]
rm(coverage, cover_nox)

cover_tot_mean <- rowMeans(cover_noxy) 
log_cover_mean <- log10(cover_tot_mean)

par(mfrow=c(1,2))
hist(log_coverage_mean, breaks = 200, xlim = c(0,5.5), main="With XY chr")
abline(v=quantile(log_coverage_mean, probs = c(seq(0.1, 1.0, by= 0.1))), col="red")
abline(v=quantile(log_coverage_mean, probs = c(0.05, 0.95), col="cyan"), lty=2)
hist(log_cover_mean, breaks = 200, xlim = c(0,5.5), main = "WhithOUT XY chr")
abline(v=quantile(log_cover_mean, probs = c(seq(0.1, 1.0, by= 0.1))), col="red")
abline(v=quantile(log_cover_mean, probs = c(0.05, 0.95), col="cyan"), lty=2)
#Looks fairly similar to the previous one 


################################  UPPER threshold  ########################################

# In order to determine where to set a threshold, we create a plot that tells us
# how many values we lose depending on where we set the threshold

coverage_histo <- unlist(cover_noxy) #This allows to inspect the whole dataframe freely
log_cover_histo <- unlist(log10(cover_noxy))
quantile(coverage_histo, probs=c(.75, .95, 1))


## UPPER THRESHOLD -- 95% quantile (5,888) to 100% (229,153)

NAs_up <- c()
thresholds_up <- c()
for(i in 0:1145){
  threshold <- 229153-i*200
  thresholds_up <- append(thresholds_up, threshold)
  Unreliable <- coverage_histo[which(coverage_histo > threshold)]
  N_Unreliable <- length(Unreliable)
  NAs_up <- append(NAs_up, N_Unreliable)
}

length(NAs_up)
length(thresholds_up) #Checking if both have the same length

par(mfrow=c(1,1))
plot(x=log10(thresholds_up), y=NAs_up, type= "l", main = "Unreliable Data by threshold", 
     xlab= "Threshold position (log)", ylab = "Unreliable Data")
abline(v=quantile(log_cover_histo, probs = c(seq(0.95, 1.0, by= 0.01))), col="red")
abline(v=quantile(log_cover_histo, probs = 0.995), col="orange")
abline(v=quantile(log_cover_histo, probs = c(seq(0, 1.0, by= 0.1))), col="green")

#Plot is far too compressed on the left side. Let's try making the x axis log10 to stretch the right part

par(mfrow=c(1,1))
plot(x=log10(thresholds_up), y=NAs_up, type= "l", main = "Unreliable Data by threshold", xlab= "threshold position", ylab = "Unreliable Data")
abline(v=quantile(log_coverage_mean, probs = c(seq(0.95, 1.0, by= 0.01))), col="green")
abline(v=quantile(log_coverage_mean, probs = 0.995), col="orange")
abline(v=quantile(coverage_tot_mean, probs = c(seq(0, 1.0, by= 0.1))), col="green")

# Let's assume the 95% quantile as our upper threshold

quantile(coverage_histo, probs=.95) # RESULT = 5.888

################################  LOWER threshold  ########################################


NAs_low <- c()
thresholds_low <- c()
for(i in 1:35){
  threshold <- i
  thresholds_low <- append(thresholds_low, threshold)
  Unreliable <- coverage_histo[which(coverage_histo < threshold)]
  N_Unreliable <- length(Unreliable)
  NAs_low <- append(NAs_low, N_Unreliable)
}


length(NAs_low)
length(thresholds_low)

par(mfrow=c(1,1))
plot( x=thresholds_low, y=NAs_low, type= "l", main = "Unreliable Data by threshold", 
      xlab= "threshold position", ylab = "Unreliable Data")
abline(v=quantile( coverage_histo, probs = c(.01, .02, .03, .04) ), 
       col=c("blue","red", "orange", "green") )

# Distribution doesn't show any "kink"; Literature recommends a minimum coverage of 30x
# (Ziller, M. J., Hansen, K. D., Meissner, A., & Aryee, M. J. (2015). 
# Coverage recommendations for methylation analysis by whole-genome bisulfite sequencing. 
# Nature methods, 12(3), 230.)

################################  APPLYING THRESHOLDS TO BETA VALUES  ########################################
# Now that we've decided where to set our coverage thresholds, every beta value whose
# coverage is either too low or too high shall be converted to NA

lower_threshold <- 30
upper_threshold <- 5888

promoters_qc <- promoters_only
for(j in 19:36){
  for(i in 1:nrow(promoters_only)) {
    if(promoters_qc [i,j]  < lower_threshold | promoters_qc [i,j] > upper_threshold){
      promoters_qc [i,j-18] <- NA
    }
  }
}
View(promoters_qc)

lower_threshold <- 30
upper_threshold <- 5888

promoters_qc <- promoters_only
for(j in 19:36){
  for(i in 1:nrow(promoters_only)) {
    if(promoters_qc [i,j]  < lower_threshold | promoters_qc [i,j] > upper_threshold){
      promoters_qc [i,j-18] <- NA
    }
  }
}

promoters_qc_beta <- promoters_qc[,c(1:18)]


###### QUALITY CONTROL OF NAs
# Let's try a similar approach as before and try to visualize the amount of NAs that there are
# in the updated beta values data frame
#     It would make sense to visualize the amount of NAs within each cohort

beta_can <- promoters_qc_beta[,c(1:9)]
beta_con <- promoters_qc_beta[,c(10:18)]

NAs_can <- rowSums(is.na(beta_can))
NAs_con <- rowSums(is.na(beta_con))
NAs_both <- rowSums(is.na(promoters_qc_beta))

NAs_all <- cbind(NAs_can, NAs_con, NAs_both)
View(NAs_all)

# We apply a similiar approach as with the coverage. We try a set of thresholds to our data
# and visualize how big the loss of data would be. 

#Cancer patients
beta_can_with_NA <- cbind(beta_can, NAs_can)

genes_left <- c()
thresholds <- c()
for(i in 0:5){
  threshold <- i
  thresholds <- append(thresholds, threshold)
  NAs_left <- nrow( NAs_all[ -which (NAs_can > threshold | NAs_con > threshold), ] )
  genes_left <- append(genes_left, NAs_left)
}

par(mfrow=c(1,1))
plot(x=thresholds, y=genes_left)

beta_all_with_NA <- cbind(beta_can, NAs_can, beta_con, NAs_con)

SEM1 <- function(x) sqrt(var(x)/length(x))

SEM2 <- function(x, na.rm=T) {
  if (na.rm) x <- na.omit(x)
  sqrt(var(x)/length(x))
}

### Quality Control part 2

data <- readRDS(file="AML_gran_list.RDS.gz")
promoters <- data$promoters[,-c(1:10)]
 
lower_threshold <- 30
upper_threshold <- 5888

# Setting unreliable data through previously calculated thresholds to NA
promoters_qc <- promoters
for(j in 19:36){
  for(i in 1:nrow(promoters)) {
    if(promoters_qc [i,j]  < lower_threshold | promoters_qc [i,j] > upper_threshold){
      promoters_qc [i,j-18] <- NA
    }
  }
}
View(promoters_qc)
# Removing genes that have more than 2 NAs in individual cohorts
beta_can <- promoters_qc[,c(1:9)]
beta_con <- promoters_qc[,c(10:18)]

NAs_can <- rowSums(is.na(beta_can))
NAs_con <- rowSums(is.na(beta_con))

Can <- cbind(beta_can, NAs_can)
Con <- cbind(beta_con, NAs_con)

beta <- cbind(Can,Con)
beta_partly_cleaned <- beta[-which(beta[,10]> 2 |beta[,20] >2),]

# Replacing NAs with the mean value of the cohort
beta_cleaned <- beta_partly_cleaned
for (i in 1:nrow(beta_cleaned)) {
   n <- beta_cleaned[i,10]
  
  if (n==1 | n==2) {
    for (j in 1:9){
      if (is.na(beta_cleaned[i,j])){
        beta_cleaned[i,j] <- rowMeans(beta_cleaned[i,1:9], na.rm = T)
      }
    }
  }
}
beta_total_cleaned <- beta_cleaned
for (i in 1:nrow(beta_total_cleaned)) {
  n <- beta_total_cleaned[i,20]
  
  if (n==1 | n==2) {
    for (j in 11:19){
      if (is.na(beta_total_cleaned[i,j])){
        beta_total_cleaned[i,j] <- rowMeans(beta_total_cleaned[i,11:19], na.rm = T)
      }
    }
  }
}
beta_ <- beta_total_cleaned[,-c(10,20)]

### Data Normalization 
promotors_normalized <- beta_
for (j in 1:18) {
  for (i in 1:nrow(promotors_normalized)) {
    promotors_normalized[[i,j]] <- log2(beta_[i,j] / (1-beta_[i,j]) )
    
  }
}
View(promotors_normalized)







    
 
