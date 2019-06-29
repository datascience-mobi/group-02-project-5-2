#Code 1: Bioinfo-05-02-QC_promoter.R

data <- readRDS(file="AML_gran_list.RDS")
annotation <- read.csv(file="sample_annotation.csv")

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




###################################### PCA #####################################

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
  paste("CON", 1:9, sep="")
)
 
View(promotors_matrix_noInf)

#The PCA is run

promotors_pca <- prcomp(t(promotors_matrix_noInf), scale = TRUE)

#promotors_pca contains three informations: x, sdev and rotation

#x contains the principal componetnts
#sdev contains how much standard deviation each component accounts for
#rotation shows how much influence(loading scores) each gene has on each PC. Positive influences push the points towards positive values and
#negative influences push the data points toward negative values.

View(promotors_pca$x)
View(promotors_pca$sdev)
View(promotors_pca$rotation)

#sdev^2 can be used to calculate how much variance each component accounts for
#We can transform that value into a percentage to better compare the PCs

pca.var <- promotors_pca$sdev^2
pca.var.per <- round(pca.var/sum(pca.var)*100, 1)

plot(pca.var.per, x = c(1 :length(pca.var.per)), type = "b", main = "Scree Plot",  xlab= "Principal Component", ylab= "Percentage of variation" )

#The first three PCs explain most of the variation within the data


#Let's program some functions to plot the samples along the PCs

ggplot.pca.cellTypeShort <- function(i, j, npar=TRUE,print=TRUE){ 
  
pca.data <- data.frame(Sample=rownames(promotors_pca$x), X=promotors_pca$x[,i], Y=promotors_pca$x[,j], PatientInfo$cellTypeShort )

ggplot(data=pca.data, aes(x=X, y=Y, label="", group = PatientInfo.cellTypeShort)) +
  geom_text() +
  xlab(paste("PC_", i," ", pca.var.per[i], "%", sep =""))  +
  ylab(paste("PC_", j," ", pca.var.per[j], "%", sep =""))  +
  theme_bw()  + 
  ggtitle("My PCA Graph") + 
 geom_point(aes(shape=PatientInfo.cellTypeShort, color=PatientInfo.cellTypeShort))

}

ggplot.pca.cellTypeGroup <- function(i, j, npar=TRUE,print=TRUE){ 
  
  pca.data <- data.frame(Sample=rownames(promotors_pca$x), X=promotors_pca$x[,i], Y=promotors_pca$x[,j], PatientInfo$cellTypeGroup )
  
  ggplot(data=pca.data, aes(x=X, y=Y, label="", group = PatientInfo.cellTypeGroup)) +
    geom_text() +
    xlab(paste("PC_", i," ", pca.var.per[i], "%", sep =""))  +
    ylab(paste("PC_", j," ", pca.var.per[j], "%", sep =""))  +
    theme_bw()  + 
    ggtitle("My PCA Graph") + 
    geom_point(aes(shape=PatientInfo.cellTypeGroup, color=PatientInfo.cellTypeGroup))
  
}

ggplot.pca.FIRST_SUBMISSION_DATE <- function(i, j, npar=TRUE,print=TRUE){ 
  
  pca.data <- data.frame(Sample=rownames(promotors_pca$x), X=promotors_pca$x[,i], Y=promotors_pca$x[,j], PatientInfo$FIRST_SUBMISSION_DATE )
  
  ggplot(data=pca.data, aes(x=X, y=Y, label="", group = PatientInfo.FIRST_SUBMISSION_DATE)) +
    geom_text() +
    xlab(paste("PC_", i," ", pca.var.per[i], "%", sep =""))  +
    ylab(paste("PC_", j," ", pca.var.per[j], "%", sep =""))  +
    theme_bw()  + 
    ggtitle("My PCA Graph") + 
    geom_point(aes(shape=PatientInfo.FIRST_SUBMISSION_DATE, color=PatientInfo.FIRST_SUBMISSION_DATE))
  
}

ggplot.pca.BIOMATERIAL_PROVIDER <- function(i, j, npar=TRUE,print=TRUE){ 
  
  pca.data <- data.frame(Sample=rownames(promotors_pca$x), X=promotors_pca$x[,i], Y=promotors_pca$x[,j], PatientInfo$BIOMATERIAL_PROVIDER )
  
  ggplot(data=pca.data, aes(x=X, y=Y, label="", group = PatientInfo.BIOMATERIAL_PROVIDER)) +
    geom_text() +
    xlab(paste("PC_", i," ", pca.var.per[i], "%", sep =""))  +
    ylab(paste("PC_", j," ", pca.var.per[j], "%", sep =""))  +
    theme_bw()  + 
    ggtitle("My PCA Graph") + 
    geom_point(aes(shape=PatientInfo.BIOMATERIAL_PROVIDER, color=PatientInfo.BIOMATERIAL_PROVIDER))
  
}

ggplot.pca.DONOR_SEX <- function(i, j, npar=TRUE,print=TRUE){ 
  
  pca.data <- data.frame(Sample=rownames(promotors_pca$x), X=promotors_pca$x[,i], Y=promotors_pca$x[,j], PatientInfo$DONOR_SEX )
  
  ggplot(data=pca.data, aes(x=X, y=Y, label="", group = PatientInfo.DONOR_SEX)) +
    geom_text() +
    xlab(paste("PC_", i," ", pca.var.per[i], "%", sep =""))  +
    ylab(paste("PC_", j," ", pca.var.per[j], "%", sep =""))  +
    theme_bw()  + 
    ggtitle("My PCA Graph") + 
    geom_point(aes(shape=PatientInfo.DONOR_SEX, color=PatientInfo.DONOR_SEX))
  
}

#Now let's plot

ggplot.pca.cellTypeShort(1,2)

ggplot.pca.cellTypeGroup(1,2)

ggplot.pca.FIRST_SUBMISSION_DATE(1,2)

ggplot.pca.BIOMATERIAL_PROVIDER(1,2)

ggplot.pca.DONOR_SEX(1,2)


#PC1 separates the control group from the AML patients which all have values below the control group ones, 
#while along PC2 the 2 AML patients' values are strogly below the control groups' values. separating them into their own group


ggplot.pca.cellTypeShort(1,3)

ggplot.pca.cellTypeGroup(1,3)

ggplot.pca.FIRST_SUBMISSION_DATE(1,3)

ggplot.pca.BIOMATERIAL_PROVIDER(1,3)

ggplot.pca.DONOR_SEX(1,3)

#PC3 strongly separates AML patients at two extremes while still preserving the control group. The two groups are characterised 
#by two different Biomaterial providers, except for one patient.

############# Logistic regression ##################

Health_Vector <- as.factor(sample_annotation$cellTypeGroup[-3])

Mvalues <- as.data.frame(t(GOI_FeatureSelection)[,1:10000])

Log_Regr_Data <- cbind(Mvalues,Health_Vector )


model <- glm(formula = Health_Vector~.,family = binomial(link = "logit"), data = Log_Regr_Data)
predict <- predict(model, newdata = Log_Regr_Data, type = "response")

#cross validation 

control = which(Log_Regr_Data$Health_Vector == "gran")

control_train = sample(control,floor(0.85*length(control)))

cancer = which(Log_Regr_Data$Health_Vector == "cancer")

cancer_train = sample(cancer,floor(0.85*length(cancer)))

train = c(control_train, cancer_train)

train_set <- Log_Regr_Data[train,] 

test_set <- Log_Regr_Data[-train,] 

model_2 <- glm(formula = Health_Vector~.,family = binomial(link = "logit"), data = train_set)
predict_2 <- predict(model_2, newdata =test_set , type = "response")


View(predict_2)  



