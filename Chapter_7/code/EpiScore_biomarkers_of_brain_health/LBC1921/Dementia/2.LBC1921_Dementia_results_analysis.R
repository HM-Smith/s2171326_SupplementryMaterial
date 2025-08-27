library(readr)
library(tidyverse)

#read in results files
Dementia_basic <- read.csv("LBC1921_dementia_84_episcores_logistic_basic_model_04072023.csv")
Dementia_full <- read.csv("LBC1921_dementia_84_episcores_logistic_full_model_no_yrsedu_04072023.csv")

Dementia_basic$SeqId <- gsub("X", "", as.character(Dementia_basic$SeqId))
Dementia_basic$SeqId <- gsub("\\.", "-", as.character(Dementia_basic$SeqId))

Dementia_full$SeqId <- gsub("X", "", as.character(Dementia_full$SeqId))
Dementia_full$SeqId <- gsub("\\.", "-", as.character(Dementia_full$SeqId))

#read in episcore ids
library(readxl)
Epi_ID <-read_excel("Index_109.xlsx")
Epi_ID <- Epi_ID %>% rename(SeqId = `Identifier (SOMAScan SeqId or Olink name)`,
                            Gene_Name = 'Gene Name')

#join results and episcore ids
Dementia_basic <- Dementia_basic %>% left_join(Epi_ID, by = "SeqId")
Dementia_full <- Dementia_full %>% left_join(Epi_ID, by = "SeqId")


#create FDR corrected column 
Dementia_basic$FDR <- p.adjust(Dementia_basic$P, method = 'BH')
Dementia_full$FDR <- p.adjust(Dementia_full$P, method = 'BH')

#filter to FDR significant results 
Dementia_basic_sig <- Dementia_basic %>% filter(FDR < 0.05)
Dementia_full_sig <- Dementia_full %>% filter(FDR < 0.05)

#check number of FDR significant associations 
dim(Dementia_basic_sig)
#[1]  3 15
dim(Dementia_full_sig)
#[1]  0 15

#filter to nominally signficant results
Dementia_basic_NOM_sig <- Dementia_basic %>% filter(P< 0.05)
Dementia_full_NOM_sig <- Dementia_full %>% filter(P < 0.05)

#check number of nominally significant results
dim(Dementia_basic_NOM_sig)
#[1]  8 15
dim(Dementia_full_NOM_sig)
#[1]  5 15


#save out files
write.csv(Dementia_basic, file = "LBC1921_dementia_84_episcores_with_ids_logistic_basic_model_05072023.csv")
write.csv(Dementia_full, file = "LBC1921_dementia_84_episcores_with_ids_logistic_full_model_no_yrsedu_05072023.csv")

write.csv(Dementia_basic_NOM_sig, file = "LBC1921_dementia_84_episcores_nominally_sig_logistic_basic_model_05072023.csv")
write.csv(Dementia_full_NOM_sig, file = "LBC1921_dementia_84_episcores_nominally_sig_logistic_full_model_no_yrsedu_05072023.csv")

write.csv(Dementia_basic_sig, file = "LBC1921_dementia_84_episcores_FDR_sig_logistic_basic_model_05072023.csv")

