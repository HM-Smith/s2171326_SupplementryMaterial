library(dplyr)

library(data.table)

#set working directory
setwd("")

#make a list of files to combine
files_intercept <- list.files(pattern = "*_G_intercept_No_Domain_LBC1936_assocs_lipid_EpiScores_basic.csv")

#bind files together
combined_intercept <- bind_rows(lapply(files_intercept, fread))

#save out combined files
write.csv(combined_intercept, file = "G_intercept_No_Domain_LBC1936_assocs_lipid_EpiScores_basic.csv")                         


#set working directory 
setwd("")

#make a list of files to combine
files_slope <- list.files(pattern = "*_G_slope_No_Domain_LBC1936_assocs_lipid_EpiScores_basic.csv")

#bind files together
combined_slope <- bind_rows(lapply(files_slope, fread))

#save out combined file
write.csv(combined_slope, file = "G_slope_No_Domain_LBC1936_assocs_lipid_EpiScores_basic.csv")

################################################################################################################################################################################################################################
################################################################################################################################################################################################################################
################################################################################################################################################################################################################################
library(dplyr)

library(data.table)

#set working directory
setwd("")

#make a list of files to combine
files_intercept <- list.files(pattern = "*_G_intercept_No_Domain_LBC1936_assocs_lipid_EpiScores_fullmodel_no_yrsedu.csv")

#bind files together 
combined_intercept <- bind_rows(lapply(files_intercept, fread))

#save out combined file
write.csv(combined_intercept, file = "G_intercept_No_Domain_LBC1936_assocs_lipid_EpiScores_fullmodel_no_yrsedu.csv")                         


#set working directory 
setwd("")

#make a list of files to combine
files_slope <- list.files(pattern = "*_G_slope_No_Domain_LBC1936_assocs_lipid_EpiScores_fullmodel_no_yrsedu.csv")

#bind files together
combined_slope <- bind_rows(lapply(files_slope, fread))

#save out combined files
write.csv(combined_slope, file = "files/G_slope_No_Domain_LBC1936_assocs_lipid_EpiScores_fullmodel_no_yrsedu_03022024.csv")


