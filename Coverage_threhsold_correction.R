data <- readRDS(file = "AML_gran_list.RDS")
promoters <- data$promoters

promoters_only <- promoters[, -c(1:10)]

beta <- promoters_only[, c(1:18)]
coverage <- promoters_only[, c(19:36)]

# Chromosomes X & Y have to be removed from the analysis
cover_nox <- coverage[-which(promoters$Chromosome == "chrX"), ]
cover_noxy <- cover_nox[-which(promoters$Chromosome == "chrY"), ]
rm(coverage, cover_nox)
beta_nox <- beta[-which(promoters$Chromosome == "chrX"), ]
beta_noxy <- beta_nox[-which(promoters$Chromosome == "chrY"), ]
rm(beta, beta_nox)

print( (nrow(beta_noxy)*100) / nrow(promoters) )


coverage_tot_mean <- rowMeans(cover_noxy) 
log_coverage_mean <- log10(cover_tot_mean)

################################  UPPER threshold  ########################################

# In order to determine where to set a threshold, we create a plot that tells us
# how many values we lose depending on where we set the threshold

coverage_histo <- unlist(cover_noxy) #This allows to inspect the whole dataframe freely
log_cover_histo <- unlist(log10(cover_noxy))

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

#How much information (= rows = promoters) did we lose at the end?
promoters_norm <- read.csv(file="promotors_normalized.csv")
(nrow(promoters_only) - nrow(promoters_norm) ) *100 / nrow(promoters_only)
#RESULT = 14.99%


# PROBLEM 1: we are not taking into account the amount of NAs that were already in the Beta values as NaN

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

# On a second thought, it would be more meaningful to check how many GENES
# we lose depending on where we set the threshold

############### LOWER THRESHOLD 

# (Ziller, M. J., Hansen, K. D., Meissner, A., & Aryee, M. J. (2015)).
# Coverage recommendations for methylation analysis by whole-genome bisulfite sequencing.
# Nature methods, 12(3), 230.)
 

lower_threshold <- 15

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

#USING FOR LOOPS

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