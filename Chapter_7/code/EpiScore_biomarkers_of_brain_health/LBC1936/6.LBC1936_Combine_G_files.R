library(dplyr)

library(data.table)

#read in and combine G intercept files then save out - basic model
setwd("")

files_intercept <- list.files(pattern = "*_G_intercept_No_Domain_LBC1936_assocs_84_EpiScores_basic_17072023.csv")

combined_intercept <- bind_rows(lapply(files_intercept, fread))

write.csv(combined_intercept, file = "G_intercept_No_Domain_LBC1936_assocs_84_EpiScores_basic_17072023.csv")                         

#read in and combine G slope files then save out - basic model
setwd("")

files_slope <- list.files(pattern = "*_G_slope_No_Domain_LBC1936_assocs_84_EpiScores_basic_17072023.csv")

combined_slope <- bind_rows(lapply(files_slope, fread))

write.csv(combined_slope, file = "G_slope_No_Domain_LBC1936_assocs_84_EpiScores_basic_17072023.csv")



#read in and combine G intercept files then save out - full model 
setwd("")

files_intercept <- list.files(pattern = "*_G_intercept_No_Domain_LBC1936_assocs_84_EpiScores_fullmodel_no_yrsedu_27102023.csv")

combined_intercept <- bind_rows(lapply(files_intercept, fread))

write.csv(combined_intercept, file = "/Cluster_Filespace/Marioni_Group/Hannah/Biomarkers_brain_health_LBC_GS/Results/LBC1936/Cognitive/combined_files/G_intercept_No_Domain_LBC1936_assocs_84_EpiScores_fullmodel_no_yrsedu_27102023.csv")                         


#read in and combine G slope files then save out -full model
setwd("")

files_slope <- list.files(pattern = "*_G_slope_No_Domain_LBC1936_assocs_84_EpiScores_fullmodel_no_yrsedu_27102023.csv")

combined_slope <- bind_rows(lapply(files_slope, fread))

write.csv(combined_slope, file = "G_slope_No_Domain_LBC1936_assocs_84_EpiScores_fullmodel_no_yrsedu_27102023.csv")


