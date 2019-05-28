#######INPUT OF THE DATA AND QC OF COVERAGE 
data <- readRDS(file="AML_gran_list.RDS")
promoters <- data$promoters;
promoters_only <- promoters[,-c(1,2,3,4,5,6,7,8,9,10)]
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

NAs_con_stat <- c(mean(NAs_con), sd(NAs_con))
NAs_con_stat
NAs_can_stat <- c(mean(NAs_can), sd(NAs_can))
NAs_can_stat
NAs_both_stat <- c(mean(NAs_both), sd(NAs_both))
NAs_both_stat
