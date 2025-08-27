library(readr)
library(tidyverse)

#read in results files
Dementia_basic <- read.csv("CompRisk_Dementia_GenScot_fineandgray_84_episcores_basic_15082023.csv")
Dementia_full <- read.csv("CompRisk_Dementia_GenScot_fineandgray_84_episcores_fullmodel_no_yrsedu_15082023.csv")

#edit episcore names so results can be combined with episcore id file
Dementia_basic$SeqId <- gsub("X", "", as.character(Dementia_basic$SeqId))
Dementia_basic$SeqId <- gsub("\\.", "-", as.character(Dementia_basic$SeqId))

Dementia_full$SeqId <- gsub("X", "", as.character(Dementia_full$SeqId))
Dementia_full$SeqId <- gsub("\\.", "-", as.character(Dementia_full$SeqId))

#read in episcore id file
library(readxl)
Epi_ID <-read_excel("Index_109.xlsx")
Epi_ID <- Epi_ID %>% rename(SeqId = `Identifier (SOMAScan SeqId or Olink name)`,
                            Gene_Name = 'Gene Name')

#combine episcore id and results files
Dementia_basic <- Dementia_basic %>% left_join(Epi_ID, by = "SeqId")
Dementia_full <- Dementia_full %>% left_join(Epi_ID, by = "SeqId")


#creat FDR adjusted column 
Dementia_basic$FDR <- p.adjust(Dementia_basic$P, method = 'BH')
Dementia_full$FDR <- p.adjust(Dementia_full$P, method = 'BH')

#filter to FDR significant results
Dementia_basic_sig <- Dementia_basic %>% filter(FDR < 0.05)
Dementia_full_sig <- Dementia_full %>% filter(FDR < 0.05)

#check number of significant associations
dim(Dementia_basic_sig)
#[1]  0 14
dim(Dementia_full_sig)
#[1]  0 14

#filter to nominally significant results
Dementia_basic_NOM_sig <- Dementia_basic %>% filter(P< 0.05)
Dementia_full_NOM_sig <- Dementia_full %>% filter(P < 0.05)

#check number of nominally significant associations 
dim(Dementia_basic_NOM_sig)
#[1]  7 14

dim(Dementia_full_NOM_sig)
#[1]  7 14

#save out results files
write.csv(Dementia_basic, file = "CompRisk_Dementia_GenScot_fineandgray_84_episcores_with_ids_basic_16082023.csv")
write.csv(Dementia_full, file = "CompRisk_Dementia_GenScot_fineandgray_84_episcores_with_ids_fullmodel_no_yrsedu_16082023.csv")

write.csv(Dementia_basic_NOM_sig, file = "CompRisk_Dementia_GenScot_fineandgray_84_episcores_nominally_sig_basic_16082023.csv")
write.csv(Dementia_full_NOM_sig, file = "CompRisk_Dementia_GenScot_fineandgray_84_episcores_nominally_sig_fullmodel_no_yrsedu_16082023.csv")



