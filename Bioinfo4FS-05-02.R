# Import dataset 
data <- readRDS(file = "C:/Users/alvar/Downloads/AML_gran_list.RDS")
annotation <- read.csv(file ="C:/Users/alvar/Downloads/sample_annotation.csv")

# Separate the dataframes into specific subgroups.
#   1st get the beta and the coverage values of the genes in 2 separate matrix
genes_beta <- data[["genes"]][,c(11:28)]
genes_coverage <- data[["genes"]][,c(29:46)]


