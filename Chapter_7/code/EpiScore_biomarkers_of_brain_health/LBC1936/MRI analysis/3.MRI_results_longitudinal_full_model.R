library(readr)
library(tidyverse)
library(dplyr)

#read in results files 
Brain <- read.csv("Totalbrain_slope_LBC1936_assocs_84_EpiScores_fullmodel_no_yrsedu_27102023.csv")
GM <- read.csv("GreyMatter_slope_LBC1936_assocs_84_EpiScores_fullmodel_no_yrsedu_27102023.csv")
NAWM <- read.csv("NAWM_slope_LBC1936_assocs_84_EpiScores_fullmodel_no_yrsedu_27102023.csv")
WMH <- read.csv("WMH_slope_LBC1936_assocs_84_EpiScores_fullmodel_no_yrsedu_27102023.csv")

#edit episcore names for merging with episcore ids 
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

#merge with episcore ids 
Brain <- Brain %>% left_join(Epi_ID, by = "SeqId")
GM <- GM %>% left_join(Epi_ID, by = "SeqId")
NAWM <- NAWM %>% left_join(Epi_ID, by = "SeqId")
WMH<- WMH %>% left_join(Epi_ID, by = "SeqId")

#make FDR corrected column
Brain$FDR <- p.adjust(Brain$P, method = 'BH')
GM$FDR <- p.adjust(GM$P, method = 'BH')
NAWM$FDR <- p.adjust(NAWM$P, method = 'BH')
WMH$FDR <- p.adjust(WMH$P, method = 'BH')

#make measure column 
Brain$Measure <- rep("Total Brain Volume", 84)
GM$Measure <- rep("Grey matter Volume", 84)
NAWM$Measure <- rep("Normal appearing White matter Volume", 84)
WMH$Measure <- rep("White matter hyperintensity volume", 84)

#filter to FDR significant results
Brain_sig <- Brain %>% filter(FDR < 0.05)
GM_sig <- GM %>% filter(FDR < 0.05)
NAWM_sig <- NAWM %>% filter(FDR < 0.05)
WMH_sig <- WMH %>% filter(FDR < 0.05)

#check FDR significant results 
dim(Brain_sig)
#[1]  0 20
dim(GM_sig)
#[1]  0 20
dim(NAWM_sig)
#[1]  0 20
dim(WMH_sig)
#[1]  0 20

#filter to nominally significant results
Brain_nom <- Brain %>% filter(P < 0.05)
GM_nom <- GM %>% filter(P < 0.05)
NAWM_nom <- NAWM %>% filter(P < 0.05)
WMH_nom <- WMH %>% filter(P < 0.05)


#check number nominally significant results 
dim(Brain_nom)
#[1]  8 20
dim(GM_nom)
#[1]  6 20
dim(NAWM_nom)
#[1]  3 20
dim(WMH_nom)
#[1]  5 20

#save out files 
write.csv(Brain, file = "Totalbrain_slope_LBC1936_assocs_84_EpiScores_with_ids_fullmodel_no_yrsedu_27102023.csv")
write.csv(GM, file = "GreyMatter_slope_LBC1936_assocs_84_EpiScores_with_ids_fullmodel_no_yrsedu_27102023.csv")
write.csv(NAWM, file = "NAWM_slope_LBC1936_assocs_84_EpiScores_with_ids_fullmodel_no_yrsedu_27102023.csv")
write.csv(WMH, file = "WMH_slope_LBC1936_assocs_84_EpiScores_withs_ids_fullmodel_no_yrsedu_27102023.csv")

write.csv(Brain_nom, file = "Totalbrain_slope_LBC1936_assocs_84_EpiScores_nominally_sig_fullmodel_no_yrsedu_27102023.csv")
write.csv(GM_nom, file = "GreyMatter_slope_LBC1936_assocs_84_EpiScores_nominally_sig_fullmodel_no_yrsedu_27102023.csv" )
write.csv(NAWM_nom, file = "NAWM_slope_LBC1936_assocs_84_EpiScores_nominally_sig_fullmodel_no_yrsedu_27102023.csv")
write.csv(WMH_nom, file = "WMH_slope_LBC1936_assocs_84_EpiScores_nominlly_sig_fullmodel_no_yrsedu_27102023.csv")



