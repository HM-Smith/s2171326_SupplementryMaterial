library(tidyverse)

#read in data
data <- read.csv("time_to_event_dementia_death_data_LBC1936_new_ascertainment_20072023.csv")

#check dementia cases
table(data$dementia_code)
#0   1
#664 107

#99999 = missing code 
summary(data$depind_w1)
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
#4    3331    5548    5443    6270   99999


#recode missing code to NA
data$depind_w1 <- ifelse(data$depind_w1 == 99999, NA, data$depind_w1)

#check recoding worked
summary(data$depind_w1)
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's
#4    3314    5511    4701    6270    6505       6

#rename data for script
data2  <- data


#make list of episcores
list_e <- colnames(data2)[24:107]
episcores = colnames(data2)[which(colnames(data2)%in% list_e)]

#create empty results dataframe 
results <- data.frame(SeqId = list_e, Odds_Ratio = NA, N = NA, ci.lower = NA, ci.upper = NA, P = NA, logOdds = NA, SE = NA)

#set rownames as episcores 
rownames(results) = results$SeqId # This will allow you to index with results[i,]

#for loop for looping over episcores and performing logistic regression 
for(i in episcores) { 
  data2$tmp = data2[,i]
  
  mod <- glm(data2$dementia_code ~ data2$tmp + scale(data2$age_baseline) + as.factor(data2$sex), family= "binomial")
  
  Beta <- exp(mod$coefficients['data2$tmp'])
  logOdds <- mod$coefficients['data2$tmp']
  P <- summary(mod)$coefficients['data2$tmp','Pr(>|z|)']
  CI.lower <- exp(confint(mod)['data2$tmp', 1])
  CI.upper <- exp(confint(mod)['data2$tmp', 2])
  N <- nobs(mod)
  SE <- summary(mod)$coefficients['data2$tmp',2]

  
  results[i,1] <- i
  results[i,2] <- Beta
  results[i,3] <- N
  results[i,4] <- CI.lower
  results[i,5] <- CI.upper
  results[i,6] <- P
  results[i,7] <- logOdds
  results[i,8] <- SE
  
  print(i)
}

#save out results 
write.csv(results, "logistic_LBC1936_Dementia_assocs_84_EpiScores_basic_model_20072023.csv")


####################################################################################################################################################################################################################
####################################################################################################################################################################################################################
#create list of episcores 
list_e <- colnames(data2)[24:107]
episcores = colnames(data2)[which(colnames(data2)%in% list_e)]

#create empty results dataframe 
results <- data.frame(SeqId = list_e, Odds_Ratio = NA, N = NA, ci.lower = NA, ci.upper = NA, P = NA, logOdds = NA, SE = NA)

#set episcores as rownames 
rownames(results) = results$SeqId # This will allow you to index with results[i,]

#for loop for looping over episcores and performing logistic regression 
for(i in episcores) { 
  data2$tmp = data2[,i]
  
  mod <- glm(data2$dementia_code ~ data2$tmp + scale(data2$age_baseline) + as.factor(data2$sex) + scale(data2$bmi_w1) + scale(data2$alcunitwk_w1) + scale(data2$depind_w1) + scale(data2$smokingScore), family= "binomial")
  
  Beta <- exp(mod$coefficients['data2$tmp'])
  logOdds <- mod$coefficients['data2$tmp']
  P <- summary(mod)$coefficients['data2$tmp','Pr(>|z|)']
  CI.lower <- exp(confint(mod)['data2$tmp', 1])
  CI.upper <- exp(confint(mod)['data2$tmp', 2])
  N <- nobs(mod)
  SE <- summary(mod)$coefficients['data2$tmp',2]
  
  
  results[i,1] <- i
  results[i,2] <- Beta
  results[i,3] <- N
  results[i,4] <- CI.lower
  results[i,5] <- CI.upper
  results[i,6] <- P
  results[i,7] <- logOdds
  results[i,8] <- SE
  
  print(i)
}

#save out results 
write.csv(results, "logistic_LBC1936_Dementia_assocs_84_EpiScores_full_model_no_yrsedu_27102023.csv")
