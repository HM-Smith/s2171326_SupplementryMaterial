library(cmprsk)
library(tidyverse)
library(readr)
library(dplyr)

#read in time to event data
data <- read.csv("time_to_event_dementia_death_data_LBC1936_new_ascertainment_20072023.csv")
data <- as.data.frame(data)

#set sex a factor
data$sex <- as.factor(data$sex)

#create status column with dementia = 1, dead = 2, unaffected = 0
data$status <- ifelse(data$dementia_code == 1, 1, ifelse(data$dead == 1, 2, 0))

#create time column with either time to dementia, time to death or time to censor 
data$time <- ifelse(data$dementia_code == 1, data$tte_dementia, data$tte_death)


#summary deprivation index, 99999 = missing code
summary(data$depind_w1)
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
#4    3331    5548    5443    6270   99999

#recode missing codes to NA
data$depind_w1 <- ifelse(data$depind_w1 == 99999, NA, data$depind_w1)

#check missing codes are NA
summary(data$depind_w1)
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's
#4    3314    5511    4701    6270    6505       6

#Select columns of interest
data1 <- data %>% select(lbc36no, age_baseline, sex, 24:107, time, status)

#check data
dim(data1)
#[1] 771 89

#subset to complete cases
data1 <- data1[complete.cases(data1),]

#check complete info 
dim(data1)
#[1] 636  89

#make list of episcores
list_e <- colnames(data1)[4:87]


#create empty dataframe 
results <- data.frame(SeqId = list_e, n = NA, beta = NA , SE = NA, P = NA, ci.upper = NA, ci.lower = NA)

#set episcores as rownames
rownames(results) = results$SeqId # This will allow you to index with results[i,]
episcores = colnames(data1)[which(colnames(data1)%in% list_e)]

#for loop that loops over episcores and performs competing risk analysis, tmp = episcore
for(i in episcores) { 
  data1$tmp = data1[,i]
  
  
  predictors <- model.matrix(~ scale(age_baseline) +
                               as.factor(sex) +
                               tmp,
                             data = data1)
  
  predictors <-  predictors[, -1]
  
  FG_crr <- crr(ftime = data1$time,
                fstatus = data1$status,
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
  
}

#save out results
write.csv(results, file = "compRisk_Dementia_LBC1936_84_EpiScores_basic_model_20072023.csv")

###########################################################################################################################################################################################
#select columns of interest
data2 <- data %>% select(lbc36no, 15:17,19:20, 24:111)

#check data
dim(data2)
#[1] 771  94

#subset to complete cases
data2 <- data2[complete.cases(data2),]

#check data of complete cases
dim(data2)
#[1] 629  94

#make list of episcores
list_e <- colnames(data2)[7:90]
episcores = colnames(data2)[which(colnames(data2)%in% list_e)]

#make empty dataframe
results <- data.frame(SeqId = list_e, n = NA, beta = NA , SE = NA, P = NA, ci.upper = NA, ci.lower = NA)

#set episcores as rownames
rownames(results) = results$SeqId # This will allow you to index with results[i,]


#for loop for looping over episcores and performing competing risk analysis 
for(i in episcores) { 
  data2$tmp = data2[,i]
  
  
  predictors <- model.matrix(~ scale(age_baseline) +
                               as.factor(sex) +
                               tmp +
                               scale(alcunitwk_w1) +
                               scale(depind_w1) +
                               scale(bmi_w1) +
                               scale(smokingScore),
                             data = data2)
  
  predictors <-  predictors[, -1]
  
  FG_crr <- crr(ftime = data2$time,
                fstatus = data2$status,
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

#save out results
write.csv(results, file = "compRisk_Dementia_LBC1936_84_EpiScores_full_model_no_yrsedu_27102023.csv")

