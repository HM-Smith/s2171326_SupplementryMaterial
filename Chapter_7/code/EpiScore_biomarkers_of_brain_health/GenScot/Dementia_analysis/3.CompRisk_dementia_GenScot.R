library(cmprsk)
library(tidyverse)
library(readr)

#load time to event data
time_to_event <- read.csv("time_to_event_data_dementia_death_GenScot_84_EpiScores_15082023.csv")

#convert to dataframe 
data <- as.data.frame(time_to_event)

#set sex as factor 
data$sex <- as.factor(data$sex)

#create status column with 0 = alive/non-affected, 1 = dementia, 2 = dead
data$status <- ifelse(data$dementia_event == 1, 1, ifelse(data$death_event == 1, 2, 0))

#create time column with either time to either dementia, death, censor 
data$time <- ifelse(data$dementia_event == 1, data$tte_dementia, data$tte_death)

#create list of EpiScores
list_e <- colnames(data)[9:92]

#create empty results dataframe 
results <- data.frame(SeqId = list_e, n = NA, beta = NA , SE = NA, P = NA, ci.upper = NA, ci.lower = NA)

#set row names as EpiScores
rownames(results) = results$SeqId # This will allow you to index with results[i,]

#for loop that loops over each EpiScore performing competing risk analysis 
episcores = colnames(data)[which(colnames(data)%in% list_e)]
for(i in episcores) { 
  data$tmp = data[,i]
  

predictors <- model.matrix(~ scale(Age_at_baseline) +
                             as.factor(sex) +
                             tmp,
                           data = data)

predictors <-  predictors[, -1]

FG_crr <- crr(ftime = data$time,
              fstatus = data$status,
              cov1 = predictors, 
              failcode = 1,
              cencode = 0)


df <- summary(FG_crr)
coeffs <- as.data.frame(df$coef)
conf <- as.data.frame(df$conf.int)
Beta <- coeffs["tmp",2]
Pvalue <- coeffs["tmp",5]
SE <- coeffs["tmp", 3]
CI_upper <- conf["tmp",4]
CI_lower <- conf["tmp",3]
n <- df$n


results[i,1] <- i
results[i,2] <- n
results[i,3] <- Beta
results[i,4] <- SE
results[i,5] <- Pvalue
results[i,6] <- CI_upper
results[i,7] <- CI_lower

print(i)

}

#save out results file
write.csv(results, file = "CompRisk_Dementia_GenScot_fineandgray_84_episcores_basic_15082023.csv")


###########################################################################################################################################################################

# read in covariate files 
alc <- read.csv("/Cluster_Filespace/Marioni_Group/Hannah/Biomarkers_brain_health_LBC_GS/GenScot_input_data/alcohol/2023-02-01_alc.csv")
smid <- read.csv("/Cluster_Filespace/Marioni_Group/Hannah/Biomarkers_brain_health_LBC_GS/GenScot_input_data/SMID/2023-02-13_simd.csv")
BMI <- read.csv("/Cluster_Filespace/Marioni_Group/GS/GS_dataset/clinical/body.csv")

#rename to sample name
names(alc)[1] <- "Sample_Name"
names(smid)[1] <- "Sample_Name"
names(BMI)[1] <- "Sample_Name"

#select columns of interest
alc <- alc[c(1,2)]
BMI <- BMI[c(1,4)]

#read in and combine EpiSmokEr scores
w1 <- readRDS("wave1_epismoker.rds")
w3 <- readRDS("wave3_epismoker.rds")
w4 <- readRDS("wave4_epismoker.rds")
bind <- rbind(w1, w3)
bind <- rbind(bind, w4) # 18779 individuals with DNAm epismoker calculated 
bind$Sample_Sentrix_ID <- row.names(bind)
bind <- bind[-1]

#join covariates and dementia data
data <- data %>% left_join(alc, by ="Sample_Name")
data <- data %>% left_join(smid, by ="Sample_Name")
data <- data %>% left_join(BMI, by ="Sample_Name")
data <- data %>% left_join(bind, by ="Sample_Sentrix_ID")

#check for missing values
table(is.na(data))
#FALSE   TRUE
#706690   1160

#loop over covariates and removed any data points 3.5 SD from the mean 
list <- c("units", "bmi", "rank","smokingScore")
for(i in list){ 
  
  cutoff1 = mean(data[,i], na.rm = T) + 3.5*sd(data[,i], na.rm = T)
  cutoff2 = mean(data[,i], na.rm = T) - 3.5*sd(data[,i], na.rm = T)
  
  data[,i][which(data[,i] > cutoff1 | data[,i] < cutoff2)] <- NA 
}

#check missing values
table(is.na(data))
#FALSE   TRUE
#706548   1302

#check dimensions of dataset
dim(data)
#[1] 7150   99

#subset to complete cases as this is required for competing risk model
data <- data[complete.cases(data), ]

#check dimensions of data
dim(data)
#[1] 5932   99

#check dementia and death cases
table(data$dementia_event, data$death_event)
     #0    1
#0 5288  502
#1   54   88



#make a list of EpiScores
list_e <- colnames(data)[9:92]

#create empty results dataframe
results <- data.frame(SeqId = list_e, n = NA, beta = NA , SE = NA, P = NA, ci.upper = NA, ci.lower = NA)

#set EpiScores as rownames
rownames(results) = results$SeqId # This will allow you to index with results[i,]

#for loop that loops over each EpiScore performing competing risk analysis 
episcores = colnames(data)[which(colnames(data)%in% list_e)]
for(i in episcores) { 
  data$tmp = data[,i]
  
  
  predictors <- model.matrix(~ scale(Age_at_baseline) +
                               as.factor(sex) +
                               tmp +
                               scale(units) +
                               scale(rank) +
                               scale(bmi) +
                               scale(smokingScore),
                             data = data)
  
  predictors <-  predictors[, -1]
  
  FG_crr <- crr(ftime = data$time,
                fstatus = data$status,
                cov1 = predictors, # for time interactions use argument cov2
                failcode = 1,
                cencode = 0)
  
  
  df <- summary(FG_crr)
  coeffs <- as.data.frame(df$coef)
  conf <- as.data.frame(df$conf.int)
  Beta <- coeffs["tmp",2]
  Pvalue <- coeffs["tmp",5]
  SE <- coeffs["tmp", 3]
  CI_upper <- conf["tmp",4]
  CI_lower <- conf["tmp",3]
  n <- df$n
  
  
  results[i,1] <- i
  results[i,2] <- n
  results[i,3] <- Beta
  results[i,4] <- SE
  results[i,5] <- Pvalue
  results[i,6] <- CI_upper
  results[i,7] <- CI_lower
  
  print(i)
  
}

#save out results file 
write.csv(results, file = "CompRisk_Dementia_GenScot_fineandgray_84_episcores_fullmodel_no_yrsedu_15082023.csv")

