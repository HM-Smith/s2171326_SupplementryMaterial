library(dplyr)

#read in phenotype-episcore correlation matrix
phenotypes_epi <- read.csv("correlation_table_for_Phenotypes_and_Elnet_EpiScore_no_glucose.csv", row.names = 1)

#rename episcores/phenotypes
names(phenotypes_epi) <- c("BMI EpiScore", "Body fat % EpiScore", "HDL cholesterol EpiScore", 
                           "Total cholesterol EpiScore", "WHR EpiScore", 
                           "BMI measured (kg/m2)", "Body fat measured (%)", "HDL cholesterol measured (mmol/L)", 
                           "Total cholesterol measured (mmol/L)", "WHR measured")
row.names(phenotypes_epi) <- c("BMI EpiScore", "Body fat % EpiScore", "HDL cholesterol EpiScore", "Total cholesterol EpiScore", "WHR EpiScore", "BMI measured (kg/m2)", "Body fat measured (%)", "HDL cholesterol measured (mmol/L)", "Total cholesterol measured (mmol/L)", "WHR measured")

#set as matrix for plotting
phenotypes_epi <- as.matrix(phenotypes_epi)


#plot correlation matrix
library(corrplot)
setwd("")
tiff("HELIOS_measured_trait_Elnet_EpiScores_correlation_heatmap.tiff", width = 9*300, height = 10*300, res = 300)
corrplot(phenotypes_epi, type = "lower", order = "hclust", tl.col = "black", addCoef.col = 'black', method = "color", main = "", number.digits = 2)
mtext("HELIOS", side = 3, line = 2, cex = 1.5)
dev.off()
