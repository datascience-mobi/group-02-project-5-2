#######INPUT OF THE DATA AND QC OF COVERAGE 
data <- readRDS(file="AML_gran_list.RDS")
promoters <- data$promoters;
promoters_only <- promoters[,-c(1:10)]
rm(promoters)

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


#Cancer patients
beta_can_with_NA <- cbind(beta_can, NAs_can)
dim(beta_can_with_NA)
beta_can_red <- beta_can[c(1:50),] # Use a reduced dataset to run preliminary coding test 


# First remove all genes that have more than 2 NAs
beta_qc_can <- beta_can_red
beta_qc_can <- beta_qc_can[-which(rowSums(is.na(beta_qc_can)) > 2), ]
View(beta_qc_can)

SEM1_can <- apply(beta_can, 1, SEM1)
View(SEM1_can)
SEM2_can <- apply(beta_can, 1, SEM2)
View(SEM2_can)

summary(SEM1_can)
summary(SEM2_can)

beta_qc2_can <- beta_qc_can
row_rm <- c()
for (i in 1:nrow(beta_qc_can)) { #this loop should be applied to a dataset in which all genes with > 2 NAs have been removed
  m <- rowMeans(beta_can_with_NA[i,], na.rm = T)
  sd2 <- sd(beta_can_with_NA[i,], na.rm = T)/sqrt(3)
    #SEM2(beta_can_with_NA[i,])
  n <- beta_can_with_NA[i,10]
  
  if(n==2){
    if(sd2 > 0.014){
      append(row_rm, i)
    } else {
      beta_qc_can[i, is.na(beta_qc_can[i,])] <- m
    }
  } else { 
    beta_qc_can[i, is.na(beta_qc_can[i,])] <- m
    }
}
row_rem <- as.vector(row_rm)
beta_qc2_can <- beta_qc2_can[-row_rem, ]

View(beta_qc2_can)

# New approach to deal with NAs
promoters_qc <- promoters
for(j in 19:36){
  for(i in 1:nrow(promoters)) {
    if(promoters_qc [i,j]  < lower_threshold | promoters_qc [i,j] > upper_threshold){
      promoters_qc [i,j-18] <- NA
    }
  }
}
View(promoters_qc)

beta_can <- promoters_qc[,c(1:9)]
beta_con <- promoters_qc[,c(10:18)]

NAs_can <- rowSums(is.na(beta_can))
NAs_con <- rowSums(is.na(beta_con))

Can <- cbind(beta_can, NAs_can)
Con <- cbind(beta_con, NAs_con)

Can_partly_cleaned <- Can[-which(rowSums(is.na(Can)) > 2), ]
View(Can_partly_cleaned)
Can_cleaned <- Can_partly_cleaned
for (i in 1:nrow(Can_cleaned)) {
  m <- rowMeans(Can_cleaned[i,], na.rm = T)
  n <- Can_cleaned[i,10]
  
  if (n==1) {
    for (j in 1:9){
      if (is.na(Can_cleaned[i,j])){
        Can_cleaned[i,j] <- rowMeans(Can_cleaned[i,], na.rm = T)
      }
    }
  }
}
View(Can_cleaned)
Can_withsd <- cbind(Can_cleaned,c(1:50995))
View(Can_withsd)
for (i in 1:nrow(Can_withsd)){
  Can_withsd[i,11] <- sd(Can_cleaned[i,], na.rm = T)/3
  
}

Can_cleaned2 <- Can_withsd

for (i in 1:nrow(Can_cleaned2)){
  n <- Can_cleaned2[i,10]
  
  if (n==2 && Can_cleaned2[i,11] < 0.15) {
    for (j in 1:9){
      if(is.na(Can_cleaned2[i,j])){
        Can_cleaned2[i,j] <- rowMeans(Can_cleaned2[i,], na.rm = T)
      }
    }
    
  }
}

View(Can_cleaned2)
