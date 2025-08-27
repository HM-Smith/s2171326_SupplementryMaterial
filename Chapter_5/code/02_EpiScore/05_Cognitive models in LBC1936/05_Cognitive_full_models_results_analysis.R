library(dplyr)

#read in cross-sectional results
LBC1936_G_intercept <- read.csv("G_intercept_No_Domain_LBC1936_assocs_lipid_EpiScores_fullmodel_no_yrsedu.csv")

#get FDR
LBC1936_G_intercept$FDR <- p.adjust(LBC1936_G_intercept$P, method = 'BH')

#filter to FDR < 0.05
LBC1936_G_intercept_SIG <- LBC1936_G_intercept %>% filter(FDR < 0.05)

#check number of FDR sig hits
dim(LBC1936_G_intercept_SIG)
#[1] 5 13

#check which traits are FDR sig
LBC1936_G_intercept_SIG$SeqId
#[1] "BMI_Epi"      "bmi_w1"       "Body_fat_Epi" "Glucose_Epi"  "WHR_Epi"

#filter to nominally sig results
LBC1936_G_intercept_nominally_sig <- LBC1936_G_intercept %>% filter(P < 0.05)

#check number of nominally sig results
dim(LBC1936_G_intercept_nominally_sig)
#[1] 5 13

#save out all results files - FDR sig, nominally sig, all results
write.csv(LBC1936_G_intercept, file = "G_intercept_No_Domain_LBC1936_assocs_lipid_EpiScores_with_ids_full_no_yrs_edu.csv")
write.csv(LBC1936_G_intercept_nominally_sig, file = "G_intercept_No_Domain_LBC1936_assocs_lipid_EpiScores_Nominally_sig_full_model_no_yrs_edu.csv")
write.csv(LBC1936_G_intercept_SIG, file = "G_intercept_No_Domain_LBC1936_assocs_lipid_EpiScores_FDR_sig_full_model_no_yrs_edu.csv")

#######################################################################################################################################################################################################################################################################################
#######################################################################################################################################################################################################################################################################################
#######################################################################################################################################################################################################################################################################################

#read in longitudinal results
LBC1936_G_slope <- read.csv("G_slope_No_Domain_LBC1936_assocs_lipid_EpiScores_fullmodel_no_yrsedu.csv")

#get FDR
LBC1936_G_slope$FDR <- p.adjust(LBC1936_G_slope$P, method = 'BH')

#filter to FDR < 0.05
LBC1936_G_slope_SIG <- LBC1936_G_slope %>% filter(FDR < 0.05)

#check number of FDR sig results
dim(LBC1936_G_slope_SIG)
# [1]  0 13

#filter to nominally sig hits
LBC1936_G_slope_nominally_sig <- LBC1936_G_slope %>% filter(P < 0.05)

#check number of nominally sig hits
dim(LBC1936_G_slope_nominally_sig)
#[1]  0 13

#save out all results file
write.csv(LBC1936_G_slope, file = "G_slope_No_Domain_LBC1936_assocs_lipid_EpiScores_with_ids_full_no_yrs_edu.csv")

