library(dplyr)

#read in projected episcores
BMI <- read.csv("LBC1936_with_BMI_Epi.csv")
HDL <- read.csv("LBC1936_with_HDL.csv")
Total <- read.csv("LBC1936_with_Total_cholesterol.csv")
Glucose <- read.csv("LBC1936_with_Glucose.csv")
Body_fat <- read.csv("LBC1936_with_Body_fat.csv")
WHR <- read.csv("LBC1936_with_WHR.csv")

#select id and scores columns
HDL <- HDL[,c(2,108)]
Total <- Total[,c(2,108)]
Glucose <- Glucose[,c(2,108)]
Body_fat <- Body_fat[,c(2,108)]
WHR<- WHR[,c(2,108)]

#join datasets together
data <- BMI %>% full_join(HDL, by = "lbc36no")
data <- data %>% full_join(Total, by = "lbc36no")
data <- data %>% full_join(Glucose, by = "lbc36no")
data <- data %>% full_join(Body_fat, by = "lbc36no")
data <- data %>% full_join(WHR, by = "lbc36no")

#save out combined dataset
write.csv(data, file = "LBC1936_with_all_lipid_EpiScores_02022024.csv")