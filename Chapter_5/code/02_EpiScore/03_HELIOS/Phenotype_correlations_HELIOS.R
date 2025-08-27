library(dplyr)
#read in phenotype correlation matrix
phenotypes <- read.csv("correlation_table_for_Phenotypes.csv", row.names = 1)

#rename phenotypes
names(phenotypes) <- c("BMI (kg/m2)", "Body fat (%)", "HDL cholesterol (mmol/L)", "Total cholesterol (mmol/L)", "WHR")
row.names(phenotypes)<- c("BMI (kg/m2)", "Body fat (%)", "HDL cholesterol (mmol/L)", "Total cholesterol (mmol/L)", "WHR")


#set as matrix for plotting
phenotypes <- as.matrix(phenotypes)


#plot correlation matrix
library(corrplot)
setwd("")
tiff("HELIOS_measured_trait_correlation_heatmap.tiff", width = 6*300, height = 7*300, res = 300)
corrplot(phenotypes, type = "lower", order = "hclust", tl.col = "black", addCoef.col = 'black', method = "color", main="", number.digits = 2)
mtext("HELIOS", side = 3, line = 2, cex = 1.5)

dev.off()
