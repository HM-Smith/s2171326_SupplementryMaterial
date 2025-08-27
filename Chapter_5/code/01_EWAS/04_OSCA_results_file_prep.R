library(dplyr)

## load bmi data 
bmi <- read.table("bmi_res.linear", header = TRUE)

## filter to "good" CpGs
cpgs_to_keep <- read.table("cpgs_tokeep.txt")

bmi_new <- bmi[bmi$Probe %in% cpgs_to_keep$V1, ]

## write out filtered to good CpG file 
write.csv(bmi_new, file = "bmi_EWAS_full_model_filtered_CpGs.csv")

## significant CpGs bmi 
bmi_sig <- bmi_new %>% filter(p < 3.6e-8)

## no of sig
dim(bmi_sig)
#[1] 57307     9

write.csv(bmi_sig, file = "bmi_EWAS_full_model_significant.csv")

#######################################################################################################################################################################################################################################


library(dplyr)

## load bmi data 
body <- read.table("body_fat_res.linear", header = T)

## filter to "good" CpGs
cpgs_to_keep <- read.table("cpgs_tokeep.txt")

body_new <- body[body$Probe %in% cpgs_to_keep$V1, ]

## write out filtered to good CpG file 
write.csv(body_new, file = "body_fat_EWAS_full_model_filtered_CpGs.csv")

## significant CpGs bmi 
body_sig <- body_new %>% filter(p < 3.6e-8)

dim(body_sig)
#[1] 29302     9


write.csv(body_sig, file = "body_fat_EWAS_full_model_significant.csv")



#######################################################################################################################################################################################################################################


library(dplyr)

## load Glucose data 
Glucose <- read.table("Glucose_res.linear", header = T)

## filter to "good" CpGs
cpgs_to_keep <- read.table("cpgs_tokeep.txt")

Glucose_new <- Glucose[Glucose$Probe %in% cpgs_to_keep$V1, ]

## write out filtered to good CpG file 
write.csv(Glucose_new, file = "Glucose_EWAS_full_model_filtered_CpGs.csv")

## significant CpGs bmi 
Glucose_sig <- Glucose_new %>% filter(p < 3.6e-8)

dim(Glucose_sig)
#[1] 460   9

write.csv(Glucose_sig, file = "Glucose_EWAS_full_model_significant.csv")

#######################################################################################################################################################################################################################################


library(dplyr)

## load Glucose data 
HDL <- read.table("HDL_cholesterol_res.linear", header = T)

## filter to "good" CpGs
cpgs_to_keep <- read.table("cpgs_tokeep.txt")

HDL_new <- HDL[HDL$Probe %in% cpgs_to_keep$V1, ]

## write out filtered to good CpG file 
write.csv(HDL_new, file = "HDL_cholesterol_EWAS_full_model_filtered_CpGs.csv")

## significant CpGs bmi 
HDL_sig <- HDL_new %>% filter(p < 3.6e-8)

dim(HDL_sig)
#[1]32288     9


write.csv(HDL_sig, file = "HDL_cholesterol_EWAS_full_model_significant.csv")

#######################################################################################################################################################################################################################################


library(dplyr)

## load Glucose data 
Total <- read.table("Total_cholesterol_res.linear", header = T)

## filter to "good" CpGs
cpgs_to_keep <- read.table("cpgs_tokeep.txt")

Total_new <- Total[Total$Probe %in% cpgs_to_keep$V1, ]

## write out filtered to good CpG file 
write.csv(Total_new, file = "Total_cholesterol_EWAS_full_model_filtered_CpGs.csv")

## significant CpGs bmi 
Total_sig <- Total_new %>% filter(p < 3.6e-8)

dim(Total_sig)
#[1] 1645    9

write.csv(Total_sig, file = "Total_cholesterol_EWAS_full_model_significant.csv")

#######################################################################################################################################################################################################################################


library(dplyr)

## load Glucose data 
whr <- read.table("whr_res.linear", header = T)

## filter to "good" CpGs
cpgs_to_keep <- read.table("cpgs_tokeep.txt")

whr_new <- whr[whr$Probe %in% cpgs_to_keep$V1, ]

## write out filtered to good CpG file 
write.csv(whr_new, file = "whr_EWAS_full_model_filtered_CpGs.csv")

## significant CpGs bmi 
whr_sig <- whr_new %>% filter(p < 3.6e-8)

dim(whr_sig)
#[1] 20622     9

write.csv(whr_sig, file = "whr_EWAS_full_model_significant.csv")
