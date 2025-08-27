library(tidyverse)
library(readxl)
library(dplyr)

#read in episcore ids 
Epi_ID <-read_excel("Index_109.xlsx")
Epi_ID <- Epi_ID %>% rename(SeqId = `Identifier (SOMAScan SeqId or Olink name)`)

#read in intercept results file - full model
LBC1936_G_intercept <- read.csv("G_intercept_No_Domain_LBC1936_assocs_84_EpiScores_fullmodel_no_yrsedu_27102023.csv")

#edit episcore names for merging with episcore ids
LBC1936_G_intercept$SeqId <- gsub("X", "", as.character(LBC1936_G_intercept$SeqId))
LBC1936_G_intercept$SeqId <- gsub("\\.", "-", as.character(LBC1936_G_intercept$SeqId))

#merge with episcore ids 
LBC1936_G_intercept <- LBC1936_G_intercept %>% left_join(Epi_ID, by = "SeqId")

#make an FDR corrected column
LBC1936_G_intercept$FDR <- p.adjust(LBC1936_G_intercept$P, method = 'BH')

#filter to FDR corrected results
LBC1936_G_intercept_SIG <- LBC1936_G_intercept %>% filter(FDR < 0.05)

#check number of FDR significant results
dim(LBC1936_G_intercept_SIG)
#[1] 0 18

#filter to nominally significant results
LBC1936_G_intercept_nominally_sig <- LBC1936_G_intercept %>% filter(P < 0.05)

#check number of nominally significant results
dim(LBC1936_G_intercept_nominally_sig)
#[1] 10 18

#save out results 
write.csv(LBC1936_G_intercept, file = "G_intercept_No_Domain_LBC1936_assocs_84_EpiScores_with_ids_full_model_no_yrsedu_27102023.csv")
write.csv(LBC1936_G_intercept_nominally_sig, file = "G_intercept_No_Domain_LBC1936_assocs_84_EpiScores_Nominally_sig_full_model_no_yrsedu_27102023.csv")

#read in slopes results file - full model
LBC1936_G_slope <- read.csv("G_slope_No_Domain_LBC1936_assocs_84_EpiScores_fullmodel_no_yrsedu_27102023.csv")

#edit episcore names to merge with episcore ids
LBC1936_G_slope$SeqId <- gsub("X", "", as.character(LBC1936_G_slope$SeqId))
LBC1936_G_slope$SeqId <- gsub("\\.", "-", as.character(LBC1936_G_slope$SeqId))

#merge with episcore ids
LBC1936_G_slope <- LBC1936_G_slope %>% left_join(Epi_ID, by = "SeqId")

#create FDR corrected column
LBC1936_G_slope$FDR <- p.adjust(LBC1936_G_slope$P, method = 'BH')

#filter to FDR significant results
LBC1936_G_slope_SIG <- LBC1936_G_slope %>% filter(FDR < 0.05)

#check number of FDR significant results
dim(LBC1936_G_slope_SIG)
#[1]  0 18

#filter to nominally significant results
LBC1936_G_slope_nominally_sig <- LBC1936_G_slope %>% filter(P < 0.05)

#check number of nominally significant results 
dim(LBC1936_G_slope_nominally_sig)
#[1]  1 18

#save out files
write.csv(LBC1936_G_slope, file = "G_slope_No_Domain_LBC1936_assocs_84_EpiScores_with_ids_full_model_no_yrsedu_30102023.csv")
write.csv(LBC1936_G_slope_nominally_sig, file = "G_slope_No_Domain_LBC1936_assocs_84_EpiScores_Nominally_sig_full_model_no_yrsedu_30102023.csv")

