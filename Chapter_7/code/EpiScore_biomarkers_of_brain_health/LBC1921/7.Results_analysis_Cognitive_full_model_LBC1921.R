library(tidyverse)
#read in results file
G_intercept_LBC1921 <- read.csv("G_intercept_No_Domain_LBC1921_assocs_84_EpiScores_full_model_no_yrsedu_21042023.csv")
G_slope_LBC1921<- read.csv("G_slope_No_Domain_LBC1921_assocs_84_EpiScores_full_model_no_yrsedu_21042023.csv")

#edit episcore name to merge with episcore ids 
G_intercept_LBC1921$SeqId <- gsub("X", "", as.character(G_intercept_LBC1921$SeqId))
G_intercept_LBC1921$SeqId <- gsub("\\.", "-", as.character(G_intercept_LBC1921$SeqId))

G_slope_LBC1921$SeqId <- gsub("X", "", as.character(G_slope_LBC1921$SeqId))
G_slope_LBC1921$SeqId <- gsub("\\.", "-", as.character(G_slope_LBC1921$SeqId))

#read in episcore ids
library(readxl)
Epi_ID <-read_excel("Index_109.xlsx")
Epi_ID <- Epi_ID %>% rename(SeqId = `Identifier (SOMAScan SeqId or Olink name)`,
                            Gene_Name = 'Gene Name')

#join with episcore ids 
G_intercept_LBC1921 <- G_intercept_LBC1921 %>% left_join(Epi_ID, by = "SeqId")
G_slope_LBC1921 <- G_slope_LBC1921 %>% left_join(Epi_ID, by = "SeqId")

#create FDR corrected column 
G_intercept_LBC1921$FDR <- p.adjust(G_intercept_LBC1921$P, method = 'BH')
G_slope_LBC1921$FDR <- p.adjust(G_slope_LBC1921$P, method = 'BH')

#filter to FDR significant results 
G_intercept_LBC1921_FDR_sig <- G_intercept_LBC1921 %>% filter(FDR < 0.05)
G_slope_LBC1921_FDR_sig <- G_slope_LBC1921 %>% filter(FDR < 0.05)

#check number of FDR significant results
dim(G_intercept_LBC1921_FDR_sig)
#[1]  5 18
dim(G_slope_LBC1921_FDR_sig)
#[1]  0 18

#filter to nominally significant results
G_intercept_LBC1921_nom_sig <- G_intercept_LBC1921 %>% filter(P< 0.05)
G_slope_LBC1921_nom_sig <- G_slope_LBC1921 %>% filter(P< 0.05)

#check number of nominally significant results
dim(G_intercept_LBC1921_nom_sig)
#[1] 19 18
dim(G_slope_LBC1921_nom_sig)
#[1]  3 18

#save out results files
write.csv(G_intercept_LBC1921, file = "G_intercept_No_Domain_LBC1921_assocs_84_EpiScores_with_ids_full_model_no_yrsedu_21042023.csv")
write.csv(G_slope_LBC1921, file = "G_slope_No_Domain_LBC1921_assocs_84_EpiScores_with_ids_full_model_no_yrsedu_21042023.csv")

write.csv(G_intercept_LBC1921_FDR_sig, file = "G_intercept_No_Domain_LBC1921_assocs_84_EpiScores_FDR_sig_full_model_no_yrsedu_21042023.csv")

write.csv(G_intercept_LBC1921_nom_sig, file = "G_intercept_No_Domain_LBC1921_assocs_84_EpiScores_Nominally_sig_full_model_no_yrsedu_21042023.csv")
write.csv(G_slope_LBC1921_nom_sig, file = "G_slope_No_Domain_LBC1921_assocs_84_EpiScores_Nominally_sig_full_model_no_yrsedu_21042023.csv" )