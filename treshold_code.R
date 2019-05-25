input_data <- readRDS(file="/Users/krank/Documents/Uni/4. Fachsemester/Bioinformatik/AMLvsGran/AML_gran_list.RDS.gz")
install.packages("rmarkdown")
sample_annotation <- read.csv(file="/Users/krank/Documents/Uni/4. Fachsemester/Bioinformatik/AMLvsGran/sample_annotation.csv")

# wir laden den Datenset von genes
Genes <- input_data$genes
Coverage<- GenesOnly[,c(19:35)]

# Plotten der Häufigkeiten bestimmter Coverage-Werte
Coverage_hist <- unlist(Coverage)
log_cover <- log10(Coverage_hist)
hist(log_cover, breaks = 200)

# schleife zum Herausfinden, wie viele Werte wir je nach Threshold verlieren und dem plotten 

tot_Unreliable <- c(1:50)
for(i in 1:50){
  threshold <- 1000000-i*5000
  Unreliable <- Coverage_hist[which(Coverage_hist > threshold)]
  N_Unreliable <- length(Unreliable)
  tot_Unreliable[i] <- N_Unreliable
  }

tot_Unreliable  

x_axis <- seq(750000, 995000, by = 5000)
plot( x=x_axis, y=tot_Unreliable, type= "l", main = "Unreliable Data by threshold", xlab= "Threshold position", ylab = "Unreliable Data")
abline(v=880000, col="red") # Our upper treshold could very well be 880.000 !

