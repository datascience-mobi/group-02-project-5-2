input_data <- readRDS(file="/Users/krank/Documents/Uni/4. Fachsemester/Bioinformatik/AMLvsGran/AML_gran_list.RDS.gz")
install.packages("rmarkdown")
sample_annotation <- read.csv(file="/Users/krank/Documents/Uni/4. Fachsemester/Bioinformatik/AMLvsGran/sample_annotation.csv")

#Split data in order to facilitate data cleanup

Genes <- input_data$genes;
Promoters <- input_data$promoters;
CpGislands <- input_data$cpgislands;

Patient_Info <- sample_annotation[, which(colnames(sample_annotation) %in% c("sampleName", "cellTypeShort", "cellTypeGroup", "tissueTypeShort", "EXPERIMENT_ID", "DONOR_AGE", "DONOR_HEALTH_STATUS", "DONOR_SEX", "DONOR_ETHNICITY", "DONOR_REGION_OF_RESIDENCE", "GENETIC_CHARACTERISTICS", "TREATMENT"))]
row.names(Patient_Info) = Patient_Info$sampleName
Patient_Info = Patient_Info[, -which(colnames(Patient_Info) %in% c("sampleName"))]

Healthy_Patients <- Patient_Info[which(Patient_Info$DONOR_HEALTH_STATUS == Healthy),]

GenesOnly <- Genes[,-c(1,2,3,4,5,6,7,8,9,10,11)]
BValues <- GenesOnly[,c(1:18)]
Coverage<- GenesOnly[,c(19:35)]

# na values in coverage 
length(is.na(Coverage)[is.na(Coverage) == T])

