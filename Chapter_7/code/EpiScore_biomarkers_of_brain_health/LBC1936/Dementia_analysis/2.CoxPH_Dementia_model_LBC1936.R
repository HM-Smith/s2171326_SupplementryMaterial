library(tidyverse)
library(readr)
library(survival)
library(kinship2)
library(coxme)
library(dplyr)

#read in data
data <- read.csv("time_to_event_dementia_death_data_LBC1936_new_ascertainment_20072023.csv")
data <- as.data.frame(data)

#set sex as factor
data$sex <- as.factor(data$sex)


#create status column, dementia = 1, unaffected = 0 
data$status <- ifelse(data$dementia_code == 1, 1, 0)

#create time column with time to dementia or time to death or time to censor
data$time <- ifelse(data$dementia_code == 1, data$tte_dementia, data$tte_death)

#dementia cases 
table(data$dementia_code)
#0   1
#664 107

#dementia cases by sex 
table(data$dementia_code, data$sex) # 2 = female, 1 = male
#  0   1
#0 345 319
#1  59  48

# 99999 = missing code 
summary(data$depind_w1)
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
#4    3331    5548    5443    6270   99999

#recode missing code to NA 
data$depind_w1 <- ifelse(data$depind_w1 == 99999, NA, data$depind_w1)

#check missing code recoded to NA
summary(data$depind_w1)
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's
#4    3314    5511    4701    6270    6505       6


#make list of episcores 
list_e <- colnames(data)[24:107]
episcores = colnames(data)[which(colnames(data)%in% list_e)]

#create empty dataframe for results
results <- data.frame(SeqId = list_e, Outcome = NA, Hazard_Ratio = NA,  ci.lower = NA,ci.upper = NA, P = NA, N_cases = NA, N_controls = NA, local = NA, global = NA)

#set episcores as rownames 
rownames(results) = results$SeqId # This will allow you to index with results[i,]

#for loop for looping over episcores and perforiming coxPH models, tmp = episcore 
for(i in episcores) { 
  data$tmp = data[,i]
  
  mod = coxph(Surv(data$time, data$status) ~ data$tmp + factor(data$sex) + scale(data$age_baseline))
  
  mod.sum <- summary(mod)$coefficients
  results[i,1] <- i
  results[i,2] <- as.character("Dementia")
  results[i,3:5]<-round(exp(cbind(coef(mod), confint(mod)))[1,1:3],2)
  results[i,6] <- mod.sum[1, 5]
  results[i,7] <- mod$nevent 
  results[i,8] <- mod$n - mod$nevent
  
  all <- cox.zph(mod) 
  p <- all$table[,"p"]
  local <- p[1]
  global <- p[4]
  
  results[i,9] <- local
  results[i,10] <- global
  
  print(i)
}

#save out results 
write.csv(results, file = "coxPH_Dementia_LBC1936_84_EpiScores_basic_model_exp_HR_20072023.csv")

###########################################################################################################################################################################################
#rename dataframe for script 
data2 <- data

#make list of episcores 
list_e <- colnames(data2)[24:107]
episcores = colnames(data2)[which(colnames(data2)%in% list_e)]

#make empty dataframe for results 
results <- data.frame(SeqId = list_e, Outcome = NA, Hazard_Ratio = NA,  ci.lower = NA,ci.upper = NA, P = NA, N_cases = NA, N_controls = NA, local = NA, global = NA)

#set episcores as rownames
rownames(results) = results$SeqId # This will allow you to index with results[i,]

#for loop for looping over episcores and performing coxPH models 
for(i in episcores) { 
  data2$tmp = data2[,i]
  
  mod = coxph(Surv(data2$time, data2$status) ~ scale(data2$tmp) + factor(data2$sex) + scale(data2$age_baseline) + scale(data2$alcunitwk_w1) + scale(data2$bmi_w1) + scale(data2$depind_w1) + scale(data2$smokingScore))
  
  mod.sum <- summary(mod)$coefficients
  results[i,1] <- i
  results[i,2] <- as.character("Dementia")
  results[i,3:5]<-round(exp(cbind(coef(mod), confint(mod)))[1,1:3],2)
  results[i,6] <- mod.sum[1, 5]
  results[i,7] <- mod$nevent 
  results[i,8] <- mod$n - mod$nevent
  
  all <- cox.zph(mod) 
  p <- all$table[,"p"]
  local <- p[1]
  global <- p[8]
  
  results[i,9] <- local
  results[i,10] <- global
  
  print(i)
}

#save out results 
write.csv(results, file = "coxPH_Dementia_LBC1936_84_EpiScores_exp_HR_full_model_no_yrsedu_27102023.csv")

