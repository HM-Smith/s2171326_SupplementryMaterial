library(dplyr)

#read in episcores
LBC_EpiScores <- read.csv("LBC1936_with_all_lipid_EpiScores_array_date_adj.csv")

#set episcores names as rownames
rownames(LBC_EpiScores) <- LBC_EpiScores$lbc36no

#select episcore columns
LBC_scores <- LBC_EpiScores[, 110:115]

#rename columns
names(LBC_scores) <- c("BMI", "HDL cholesterol", "Total cholesterol", "Glucose", "Body fat %", "WHR")


# create correlation matrix 
res <- cor(LBC_scores)


write.csv(res, file = "/Cluster_Filespace/Marioni_Group/Hannah/Lipid_Related_Trait_EWAS_GS20K/EpiScore_correlations/EpiScore_adj_correlations_table_07032024.csv")
 
library(corrplot)
setwd("/Cluster_Filespace/Marioni_Group/Hannah/Lipid_Related_Trait_EWAS_GS20K/EpiScore_correlations/")
tiff("EpiScore_adj_correlation_heatmap_07032024.tiff", width = 6*300, height = 7*300, res = 300)
corrplot(res, type = "upper", order = "hclust", tl.col = "black", addCoef.col = 'black', method = "color", main = "", number.digits = 2)
mtext("LBC1936", side = 3, line = 2, cex = 1.5)
dev.off()


###############################################################################################################################
data <- LBC_EpiScores[, c(7,9, 11)]

names(data) <- c("BMI (kg/m2)", "Total cholesterol (mmol/L)", "HDL cholesterol (mmol/L)")
res <- cor(data, use="pairwise.complete.obs")

#save out correlation matrix
write.csv(res, file ="/Cluster_Filespace/Marioni_Group/Hannah/Lipid_Related_Trait_EWAS_GS20K/measure_episcore_comparison_LBC1936/measured_trait_corrs_lbc1936_07032024.csv")


library(corrplot)

#plot correlation matrix
setwd("")
tiff("measured_trait_LBC1936_correlation_heatmap.tiff", width = 6*300, height = 7*300, res = 300)
corrplot(res, type = "upper", order = "hclust", tl.col = "black", addCoef.col = 'black', method = "color", main="", number.digits = 2)
mtext("LBC1936", side = 3, line = 2, cex = 1.5)
dev.off()


###########################################################################################################################
library(dplyr)

#read in episcores and measured traits
LBC_EpiScores <- read.csv("LBC1936_with_all_lipid_EpiScores_array_date_adj.csv")

#set episcore and measure trait names as rownames
rownames(LBC_EpiScores) <- LBC_EpiScores$lbc36no

#select episcore and measured trait columns
LBC_scores <- LBC_EpiScores[, c(7,9,11,110:112)]

#rename columns
names(LBC_scores) <- c("BMI measured (kg/m2)", "Total cholesterol measured (mmol/L)", "HDL cholesterol measured (mmol/L)", "BMI EpiScore", "HDL cholesterol EpiScore", "Total cholesterol EpiScore")


#create correlation matrix
res <- cor(LBC_scores, use='pairwise.complete.obs')

#save out correlation matrix 
write.csv(res, file = "EpiScore_adj_measured_trait_correlations_table.csv")

library(corrplot)
setwd("")
tiff("EpiScore_adj_measured_trait_correlation_heatmap.tiff", width = 6*300, height = 7*300, res = 300)
corrplot(res, type = "upper", order = "hclust", tl.col = "black", addCoef.col = 'black', method = "color", main = "", number.digits = 2)
mtext("LBC1936", side = 3, line = 2, cex = 1.5)

dev.off()
