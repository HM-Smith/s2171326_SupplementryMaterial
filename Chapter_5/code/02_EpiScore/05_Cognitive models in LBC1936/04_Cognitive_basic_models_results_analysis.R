library(dplyr)
#read in cross-sectional results file
LBC1936_G_intercept <- read.csv("G_intercept_No_Domain_LBC1936_assocs_lipid_EpiScores_basic.csv")

#get FDR 
LBC1936_G_intercept$FDR <- p.adjust(LBC1936_G_intercept$P, method = 'BH')

#filter results by FDR < 0.05
LBC1936_G_intercept_SIG <- LBC1936_G_intercept %>% filter(FDR < 0.05)

#check num of FDR sig results
dim(LBC1936_G_intercept_SIG)
#[1] 8 13

#check which traits were FDR significant
LBC1936_G_intercept_SIG$SeqId
#[1] "bld_choles_w1"  "bld_hdlchol_w1" "BMI_Epi"        "bmi_w1"
#[5] "Body_fat_Epi"   "Glucose_Epi"    "HDL_Epi"        "WHR_Epi"


#filter to nominally sig results (p < 0.05)
LBC1936_G_intercept_nominally_sig <- LBC1936_G_intercept %>% filter(P < 0.05)

#check number of nominally sig results
dim(LBC1936_G_intercept_nominally_sig)
#[1] 6 13

#save out FDR sig, nominally sig, and full result files
write.csv(LBC1936_G_intercept, file = "G_intercept_No_Domain_LBC1936_assocs_lipid_EpiScores_with_ids_basic.csv")
write.csv(LBC1936_G_intercept_SIG, file = "G_intercept_No_Domain_LBC1936_assocs_lipid_EpiScores_FDR_sig_basic.csv")
write.csv(LBC1936_G_intercept_nominally_sig, file = "G_intercept_No_Domain_LBC1936_assocs_lipid_EpiScores_Nominally_sig_basic.csv")

#######################################################################################################################################################################################################################################################################################
#######################################################################################################################################################################################################################################################################################
#######################################################################################################################################################################################################################################################################################

#read in longitudinal results file
LBC1936_G_slope <- read.csv("/Cluster_Filespace/Marioni_Group/Hannah/Lipid_Related_Trait_EWAS_GS20K/Cognitive_results/combined_files/G_slope_No_Domain_LBC1936_assocs_lipid_EpiScores_basic_03022024.csv")

#get FDR 
LBC1936_G_slope$FDR <- p.adjust(LBC1936_G_slope$P, method = 'BH')
#filter by FDR < 0.05
LBC1936_G_slope_SIG <- LBC1936_G_slope %>% filter(FDR < 0.05)

#check number of FDR sig results
dim(LBC1936_G_slope_SIG)
# [1]  0 13

#filter to nominally significant results
LBC1936_G_slope_nominally_sig <- LBC1936_G_slope %>% filter(P < 0.05)

# check number of nominally significant results
dim(LBC1936_G_slope_nominally_sig)
#[1]  0 13

#save out results file
write.csv(LBC1936_G_slope, file = "G_slope_No_Domain_LBC1936_assocs_lipid_EpiScores_with_ids_basic.csv")

