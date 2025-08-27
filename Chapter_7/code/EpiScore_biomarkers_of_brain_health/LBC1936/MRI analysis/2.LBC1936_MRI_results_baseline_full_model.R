library(readr)
library(tidyverse)
library(dplyr)

#read in results files 
Brain <- read.csv("Totalbrain_intercept_LBC1936_assocs_84_EpiScores_fullmodel_no_yrsedu_27102023.csv")
GM <- read.csv("GreyMatter_intercept_LBC1936_assocs_84_EpiScores_fullmodel_no_yrsedu_27102023.csv")
NAWM <- read.csv("NAWM_intercept_LBC1936_assocs_84_EpiScores_fullmodel_no_yrsedu_27102023.csv")
WMH <- read.csv("/WMH_intercept_LBC1936_assocs_84_EpiScores_fullmodel_no_yrsedu_27102023.csv")


#edit episcore names to merge with episcore ids 
Brain$SeqId <- gsub("X", "", as.character(Brain$SeqId))
Brain$SeqId <- gsub("\\.", "-", as.character(Brain$SeqId))

GM$SeqId <- gsub("X", "", as.character(GM$SeqId))
GM$SeqId <- gsub("\\.", "-", as.character(GM$SeqId))

NAWM$SeqId <- gsub("X", "", as.character(NAWM$SeqId))
NAWM$SeqId <- gsub("\\.", "-", as.character(NAWM$SeqId))

WMH$SeqId <- gsub("X", "", as.character(WMH$SeqId))
WMH$SeqId <- gsub("\\.", "-", as.character(WMH$SeqId))


#read in episcore ids 
library(readxl)
Epi_ID <-read_excel("Index_109.xlsx")
Epi_ID <- Epi_ID %>% rename(SeqId = `Identifier (SOMAScan SeqId or Olink name)`,
                            Gene_Name = 'Gene Name')

#join with episcore ids 
Brain <- Brain %>% left_join(Epi_ID, by = "SeqId")
GM <- GM %>% left_join(Epi_ID, by = "SeqId")
NAWM <- NAWM %>% left_join(Epi_ID, by = "SeqId")
WMH<- WMH %>% left_join(Epi_ID, by = "SeqId")

#create FDR corrected column 
Brain$FDR <- p.adjust(Brain$P, method = 'BH')
GM$FDR <- p.adjust(GM$P, method = 'BH')
NAWM$FDR <- p.adjust(NAWM$P, method = 'BH')
WMH$FDR <- p.adjust(WMH$P, method = 'BH')

#create measure column 
Brain$Measure <- rep("Total Brain Volume", 84)
GM$Measure <- rep("Grey matter Volume", 84)
NAWM$Measure <- rep("Normal appearing White matter Volume", 84)
WMH$Measure <- rep("White matter hyperintensity volume", 84)

#filter to FDR significant results 
Brain_sig <- Brain %>% filter(FDR < 0.05)
GM_sig <- GM %>% filter(FDR < 0.05)
NAWM_sig <- NAWM %>% filter(FDR < 0.05)
WMH_sig <- WMH %>% filter(FDR < 0.05)

#check number of FDR significant results 
dim(Brain_sig)
#[1]  0 20
dim(GM_sig)
#[1]  1 20
dim(NAWM_sig)
#[1]  0 20
dim(WMH_sig)
#[1]  0 20


#filter to nominally significant results 
Brain_nom <- Brain %>% filter(P < 0.05)
GM_nom <- GM %>% filter(P < 0.05)
NAWM_nom <- NAWM %>% filter(P < 0.05)
WMH_nom <- WMH %>% filter(P < 0.05)


#check number of nominally significant results
dim(Brain_nom)
#[1]  6 20
dim(GM_nom)
#[1]  8 20
dim(NAWM_nom)
#[1] 12 20
dim(WMH_nom)
#[1]  6 20

  
#save out results files 
write.csv(Brain, file = "Totalbrain_intercept_LBC1936_assocs_84_EpiScores_with_ids_fullmodel_no_yrsedu_27102023.csv")
write.csv(GM, file = "GreyMatter_intercept_LBC1936_assocs_84_EpiScores_with_ids_fullmodel_no_yrsedu_27102023.csv")
write.csv(NAWM, file = "NAWM_intercept_LBC1936_assocs_84_EpiScores_woth_ids_fullmodel_no_yrsedu_27102023.csv")
write.csv(WMH, file = "WMH_intercept_LBC1936_assocs_84_EpiScores_with_ids_fullmodel_no_yrsedu_27102023.csv")

write.csv(GM_sig, file = "GreyMatter_intercept_LBC1936_assocs_84_EpiScores_FDR_sig_fullmodel_no_yrsedu_27102023.csv")


write.csv(Brain_nom, file = "Totalbrain_intercept_LBC1936_assocs_84_EpiScores_Nominally_sig_fullmodel_no_yrsedu_27102023.csv")
write.csv(GM_nom, file = "GreyMatter_intercept_LBC1936_assocs_84_EpiScores_Nominally_sig_fullmodel_no_yrsedu_27102023.csv" )
write.csv(NAWM_nom, file = "NAWM_intercept_LBC1936_assocs_84_EpiScores_Nominally_sig_fullmodel_no_yrsedu_27102023.csv")
write.csv(WMH_nom, file = "WMH_intercept_LBC1936_assocs_84_EpiScores_Nominally_sig_fullmodel_no_yrsedu_27102023.csv")









GM_sig$Gene_Name = factor(GM_sig$Gene_Name, levels=unique(GM_sig$Gene_Name[order(GM_sig$beta)]))
NAWM_sig$Gene_Name = factor(NAWM_sig$Gene_Name, levels=unique(NAWM_sig$Gene_Name[order(NAWM_sig$beta)]))




library(ggplot2)
setwd("/Cluster_Filespace/Marioni_Group/Hannah/Biomarkers_brain_health_LBC_GS/Results/LBC1936/MRI/Plots/")
tiff("GM_Volume_Baseline_LBC1936_fullmodel_no_yrsedu_15022023.tiff")
ggplot(GM_sig, aes(x = beta, y = Gene_Name)) + 
  geom_point(size = 4.5, color = 'palegreen4') +
  geom_errorbarh(aes(xmax = ci.upper, xmin = ci.lower), size = .9, height = 
                   .4, color = 'palegreen4') +
  geom_vline(aes(xintercept = 0), size = .25, linetype = "dashed") + 
  theme_classic()+ 
  ylab("EpiScore") +
  xlab("Beta") +
  ggtitle("Grey matter volume at Baseline")  
dev.off()



library(ggplot2)
setwd("/Cluster_Filespace/Marioni_Group/Hannah/Biomarkers_brain_health_LBC_GS/Results/LBC1936/MRI/Plots/")
tiff("NAWM_Volume_Baseline_LBC1936_fullmodel_no_yrsedu_15022023.tiff")
ggplot(NAWM_sig, aes(x = beta, y = Gene_Name)) + 
  geom_point(size = 4.5, color = 'slateblue') +
  geom_errorbarh(aes(xmax = ci.upper, xmin = ci.lower), size = .9, height = 
                   .4, color = 'slateblue') +
  geom_vline(aes(xintercept = 0), size = .25, linetype = "dashed") + 
  theme_classic()+ 
  ylab("EpiScore") +
  xlab("Beta") +
  ggtitle("Normal appearing white matter volume at Baseline")  
dev.off()



library(cowplot)


a <- ggplot(GM_sig, aes(x = beta, y = Gene_Name)) + 
  geom_point(size = 4.5, color = 'palegreen4') +
  geom_errorbarh(aes(xmax = ci.upper, xmin = ci.lower), size = .9, height = 
                   .4, color = 'palegreen4') +
  geom_vline(aes(xintercept = 0), size = .25, linetype = "dashed") + 
  theme_classic()+ 
  ylab("EpiScore") +
  xlab("Beta") +
  ggtitle("Grey matter volume at Baseline")  

b <- ggplot(NAWM_sig, aes(x = beta, y = Gene_Name)) + 
  geom_point(size = 4.5, color = 'slateblue') +
  geom_errorbarh(aes(xmax = ci.upper, xmin = ci.lower), size = .9, height = 
                   .4, color = 'slateblue') +
  geom_vline(aes(xintercept = 0), size = .25, linetype = "dashed") + 
  theme_classic()+ 
  ylab("EpiScore") +
  xlab("Beta") +
  ggtitle("Normal appearing white matter volume at Baseline")  




setwd("/Cluster_Filespace/Marioni_Group/Hannah/Biomarkers_brain_health_LBC_GS/Results/LBC1936/MRI/Plots/")
tiff("ALL_MRI_Volume_Baseline_LBC1936_full_model_no_yrsedu_15022023.tiff", width = 900, height = 900)
plot_grid(a,b)
dev.off()





