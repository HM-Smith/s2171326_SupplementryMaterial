library(readr)
library(tidyverse)
library(dplyr)

#read in results files 
Dementia_basic <- read.csv("compRisk_Dementia_LBC1936_84_EpiScores_basic_model_20072023.csv")
Dementia_full <- read.csv("compRisk_Dementia_LBC1936_84_EpiScores_full_model_no_yrsedu_27102023.csv")


#edit episcore names for merging with episcore ids 
Dementia_basic$SeqId <- gsub("X", "", as.character(Dementia_basic$SeqId))
Dementia_basic$SeqId <- gsub("\\.", "-", as.character(Dementia_basic$SeqId))

Dementia_full$SeqId <- gsub("X", "", as.character(Dementia_full$SeqId))
Dementia_full$SeqId <- gsub("\\.", "-", as.character(Dementia_full$SeqId))

#read in episcore ids 
library(readxl)
Epi_ID <-read_excel("Index_109.xlsx")
Epi_ID <- Epi_ID %>% rename(SeqId = `Identifier (SOMAScan SeqId or Olink name)`,
                            Gene_Name = 'Gene Name')

#join with episcore ids 
Dementia_basic <- Dementia_basic %>% left_join(Epi_ID, by = "SeqId")
Dementia_full <- Dementia_full %>% left_join(Epi_ID, by = "SeqId")

#create FDR corrected column 
Dementia_basic$FDR <- p.adjust(Dementia_basic$P, method = 'BH')
Dementia_full$FDR <- p.adjust(Dementia_full$P, method = 'BH')

#filter to FDR significant results 
Dementia_basic_sig <- Dementia_basic %>% filter(FDR < 0.05)
Dementia_full_sig <- Dementia_full %>% filter(FDR < 0.05)

#check number of FDR significant results
dim(Dementia_basic_sig)
#[1]  0 14
dim(Dementia_full_sig)
#[1]  0 14

#filter to nominally significant results
Dementia_basic_NOM_sig <- Dementia_basic %>% filter(P< 0.05)
Dementia_full_NOM_sig <- Dementia_full %>% filter(P < 0.05)

#check number of nominally significant results 
dim(Dementia_basic_NOM_sig)
#[1]  5 14
dim(Dementia_full_NOM_sig)
#[1]  3 14

#save out results 
write.csv(Dementia_basic, file = "compRisk_Dementia_LBC1936_84_EpiScores_with_ids_basic_model_20072023.csv")
write.csv(Dementia_full, file = "compRisk_Dementia_LBC1936_84_EpiScores_with_ids_full_model_no_yrsedu_27102023.csv")

write.csv(Dementia_basic_NOM_sig, file = "compRisk_Dementia_LBC1936_84_EpiScores_nominally_sig_basic_model_20072023.csv")
write.csv(Dementia_full_NOM_sig, file = "compRisk_Dementia_LBC1936_84_EpiScores_nominally_sig_full_model_no_yrsedu_27102023.csv")



