library(dplyr)

#read in HELIOS EpiScore correlation matrix
EpiScores <- read.csv("Helios_correlation_table_for_six_EpiScores_wholecohort.csv")

#set episcore names as rownames
rownames(EpiScores) <- EpiScores[,1]

#remove first column 
EpiScores <- EpiScores[,-1]

#rename episcores
names(EpiScores) <- c("BMI", "Body fat %", "Glucose", "HDL cholesterol", "Total cholesterol", "WHR")

#set as matrix for plotting
EpiScores <- as.matrix(EpiScores)


#plot correlation matrix 
library(corrplot)
setwd("")
tiff("HELIOS_EpiScore_correlation_heatmap.tiff", width = 6*300, height = 7*300, res = 300)
corrplot(EpiScores, type = "lower", order = "hclust", tl.col = "black", addCoef.col = 'black', method = "color", main = "", number.digits = 2)
mtext("HELIOS", side = 3, line = 2, cex = 1.5)
dev.off()
