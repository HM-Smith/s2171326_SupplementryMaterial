library(dplyr)

#read in EWAS cat datasets and combine
studies <- read.delim("ewascatalog-studies.txt", sep = "\t")
results <- read.delim("ewascatalog-results.txt", sep = "\t")
EWAS_cat <- results %>% left_join(studies, by = "StudyID")

#make a list of trait names and save out for browsing 
Trait_names <- unique(EWAS_cat$Trait)
Trait_names <- as.data.frame(Trait_names)
write.csv(Trait_names, file = "/Trait_names.csv")



#filter trait names for BMI
BMI_ewas_cat <- EWAS_cat %>% filter(Trait == "BMI" | Trait == "Body Mass Index" | Trait == "bmi" | 
                                      Trait == "Body mass index" | Trait == "Pre-pregnancy BMI" | 
                                      Trait == "Pre-pregnancy body mass index" | Trait == "Body mass index change" |
                                      Trait == "body mass index" | Trait == "Paternal body mass index" | Trait == "maternal pre-pregnancy body mass index")
dim(BMI_ewas_cat)
#[1] 6427   32

#save out BMI studies
write.csv(BMI_ewas_cat, file = "BMI_EWAS_cat.csv") 

#filter trait names for WHR
WHR_ewas_cat <- EWAS_cat %>% filter(Trait == "waist circumference-to-hip ratio")

dim(WHR_ewas_cat)
#[1] 10 32

#save out WHR studies
write.csv(WHR_ewas_cat, file = "WHR.csv")

#filter trait names for body fat
bodyfat_ewas_Cat <- EWAS_cat %>% filter(Trait == ) # couldn't find any 

#filter trait names for glucose studies
Glucose_ewas_cat <- EWAS_cat %>% filter(Trait == "2-hour glucose" | Trait == "Fasting glucose"| 
                                          Trait == "fasting glucose" | Trait == "1-hour glucose" |
                                          Trait == "Mid-pregnancy Fasting Plasma Glucose" | 
                                          Trait == "Mid-pregnancy Fasting Plasma Glucose Post 2h 75g Oral Glucose Tolerance Test" |
                                          Trait == "Glucose" | Trait == "maternal glucose" | Trait == "oral glucose tolerance test")
#check data dimensions 
dim(Glucose_ewas_cat)
#[1] 3162   32

#save out file
write.csv(Glucose_ewas_cat, file = "Glucose_ewas_cat.csv")

#filter trait to HDL cholesterol
HDL_ewas_cat <- EWAS_cat %>% filter(Trait == "HDL cholesterol" | 
                                      Trait == "High-density lipoprotein cholesterol" | 
                                      Trait == "Serum high-densitty lipoprotein cholesterol" | 
                                      Trait == "Cholesterol esters in large HDL" |
                                      Trait == "Cholesterol esters in medium HDL" |
                                      Trait == "Cholesterol esters in small HDL" |
                                      Trait == "Cholesterol esters in very large HDL" |
                                      Trait == "Concentration of large HDL particles" |
                                      Trait == "Concentration of medium HDL particles" |
                                      Trait == "Concentration of small HDL particles" |
                                      Trait == "Concentration of very large HDL particles" |
                                      Trait == "Free cholesterol in large HDL" |
                                      Trait == "Free cholesterol in medium HDL" |
                                      Trait == "Free cholesterol in small HDL" |
                                      Trait == "Free cholesterol in very large HDL" |
                                      Trait == "Total cholesterol in HDL" |
                                      Trait == "Total cholesterol in HDL2" |
                                      Trait == "Total cholesterol in HDL3" |
                                      Trait == "Total cholesterol in large HDL")

#check data dimensions 
dim(HDL_ewas_cat)
#[1] 1026   32

#save out file
write.csv(HDL_ewas_cat, file = "HDL_ewas_cat.csv")

#filter trait to total cholesterol
Totat_chol_ewas_cat <- EWAS_cat %>% filter(Trait == "Total cholesterol" | Trait == "total cholesterol" | Trait == "Serum total cholesterol")

#check data dimensions 
dim(Totat_chol_ewas_cat)
#[1] 250  32

#save out file 
write.csv(Totat_chol_ewas_cat, file = "Total_cholesterol_EWAS_cat.csv")