library(readr)
library(tidyverse)

#read in results files
Dementia_basic <- read.csv("Logistic_Dementia_GenScot_assocs_84_EpiScores_basic_model_updated_15082023.csv")
Dementia_full <- read.csv("Logistic_Dementia_GenScot_assocs_84_EpiScores_full_model_no_yrsedu_updated_15082023.csv")

#edit episcore names to merge with episcore ids
Dementia_basic$SeqId <- gsub("X", "", as.character(Dementia_basic$SeqId))
Dementia_basic$SeqId <- gsub("\\.", "-", as.character(Dementia_basic$SeqId))

Dementia_full$SeqId <- gsub("X", "", as.character(Dementia_full$SeqId))
Dementia_full$SeqId <- gsub("\\.", "-", as.character(Dementia_full$SeqId))

#load episcore ids 
library(readxl)
Epi_ID <-read_excel("Index_109.xlsx")
Epi_ID <- Epi_ID %>% rename(SeqId = `Identifier (SOMAScan SeqId or Olink name)`,
                            Gene_Name = 'Gene Name')

#join results files and episcore ids 
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
#[1]  0 15
dim(Dementia_full_sig)
#[1]  0 15

#filter to nominally significant results
Dementia_basic_NOM_sig <- Dementia_basic %>% filter(P< 0.05)
Dementia_full_NOM_sig <- Dementia_full %>% filter(P < 0.05)

#check number of nominally significant results
dim(Dementia_basic_NOM_sig)
#[1]  9 15
dim(Dementia_full_NOM_sig)
#[1]  6 15

#save out results files
write.csv(Dementia_basic, file = "logistic_Dementia_GenScot_84_episcores_with_ids_basic_16082023.csv")
write.csv(Dementia_full, file =  "logistic_Dementia_GenScot_84_episcores_with_ids_Full_model_no_yrsedu_16082023.csv")

write.csv(Dementia_basic_NOM_sig, file = "logistic_Dementia_GenScot_84_episcores_nominally_sig_basic_16082023.csv")
write.csv(Dementia_full_NOM_sig, file = "logistic_Dementia_GenScot_84_episcores_nominally_sig_Full_model_no_yrsedu_16082023.csv" )



