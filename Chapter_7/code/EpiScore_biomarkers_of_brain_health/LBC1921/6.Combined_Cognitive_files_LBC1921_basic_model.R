library(dplyr)

library(data.table)

#combine and save out G intercept results files
setwd("")

files_intercept <- list.files(pattern ="*_G_intercept_No_Domain_LBC1921_assocs_84_EpiScores_basic_21042023.csv")

combined_intercept <- bind_rows(lapply(files_intercept, fread))

write.csv(combined_intercept, file = "G_intercept_No_Domain_LBC1921_assocs_84_EpiScores_basic_model_21042023.csv")                         


#combine and save out G slope results files
setwd("")

files_slope <- list.files(pattern = "*_G_slope_No_Domain_LBC1921_assocs_84_EpiScores_basic_21042023.csv")

combined_slope <- bind_rows(lapply(files_slope, fread))

write.csv(combined_slope, file = "G_slope_No_Domain_LBC1921_assocs_84_EpiScores_basic_model_21042023.csv")                         








