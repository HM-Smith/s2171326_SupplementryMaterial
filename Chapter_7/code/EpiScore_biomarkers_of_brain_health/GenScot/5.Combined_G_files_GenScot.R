library(dplyr)

library(data.table)

#set working directory 
setwd("")

#make list of basic model results
files_intercept <- list.files(pattern = "*_G__Generation_Scot_assocs_84_EpiScores_basic_20012023.csv")

#combine files in list
combined_intercept <- bind_rows(lapply(files_intercept, fread))

#write out combined file
write.csv(combined_intercept, file = "GS20K_G_intercept_assocs_84_EpiScores_basic_20012023.csv")


#set working directory 
setwd("")

#make list of full model results
files_intercept <- list.files(pattern = "*_G__Generation_Scot_assocs_84_EpiScores_fullmodel_no_yrsedu_13022023.csv")

#combined files in list
combined_intercept <- bind_rows(lapply(files_intercept, fread))

#write out combined file
write.csv(combined_intercept, file = "GS20K_G_intercept_assocs_84_EpiScores_fullmodel_no_yrsedu_14022023.csv")

