library(tidyverse)
library(readxl)

#read in episcore ids 
Epi_ID <-read_excel("Index_109.xlsx")
Epi_ID <- Epi_ID %>% rename(SeqId = `Identifier (SOMAScan SeqId or Olink name)`)

#read in full model results
GenScot_G_intercept <- read.csv("GS20K_G_intercept_assocs_84_EpiScores_fullmodel_no_yrsedu_14022023.csv")

#edit episcore names to join with episcore ids
GenScot_G_intercept$SeqId <- gsub("X", "", as.character(GenScot_G_intercept$SeqId))
GenScot_G_intercept$SeqId <- gsub("\\.", "-", as.character(GenScot_G_intercept$SeqId))

#join episcore ids and results
GenScot_G_intercept <- GenScot_G_intercept %>% left_join(Epi_ID, by = "SeqId")

#create FDR corrected column 
GenScot_G_intercept$FDR <- p.adjust(GenScot_G_intercept$P, method = 'BH')

#filter to FDR significant results
GenScot_G_intercept_SIG <- GenScot_G_intercept %>% filter(FDR < 0.05)

#check number of FDR signficant associations 
dim(GenScot_G_intercept_SIG)
#40 19

#filter to nominally sigificant results
GenScot_G_intercept_nominally_sig <- GenScot_G_intercept %>% filter(P < 0.05)

#check number of nominally significant associations 
dim(GenScot_G_intercept_nominally_sig)
#[1] 49 19

#save out files
write.csv(GenScot_G_intercept, file = "GS20K_G_intercept_assocs_84_EpiScores_with_ids_full_model_no_yrsedu_31032023.csv")
write.csv(GenScot_G_intercept_SIG, file = "GS20K_G_intercept_assocs_84_EpiScores_Nominally_sig_full_model_no_yrsedu_31032023.csv")
