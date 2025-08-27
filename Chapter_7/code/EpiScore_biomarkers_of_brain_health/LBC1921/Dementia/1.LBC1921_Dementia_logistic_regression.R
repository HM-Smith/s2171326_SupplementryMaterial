library(tidyverse)

#read in dementia data
data <- read.csv("LBC1921_WAVE1_cohort_data_KORA_EpiScore_Adjusted_31012023.csv")

#check dementia case numbers
table(data$dement_consensus)

      #M   N   Y
#439   7  13 110

#check dementia subtype numbers 
table(data$dement_subtype)

#          MCI    POSSIBLE MCI POSSIBLE PDD  POSSIBLE VD  PROBABLE AD
# 449      2            1            1            7           38
#PROBABLE MIX PROBABLE PSP  PROBABLE VD    UNKNOWN
#      9            1           25           36


#select columns of interest 
data2  <- data %>% select(3, 32, 41:42, 49, 102:189, 193)

#read in EpiSmokEr scores and rename id column name to basename  
epismoke <- readRDS("lbc_epismoker.rds") 
names(epismoke)[1] <- "Basename"

#combine dementia data with epismoker scores
data2 <- data2 %>% left_join(epismoke, by = "Basename")

#create status column where dementia = 1, and non-affected = 0
data2$status <- ifelse(data2$dement_consensus == "Y", 1, 0)

#check dementia cases
table(data2$dement_consensus)
#0   1
#459 110

#recode possible dementia cases tp NA
data2$status <- ifelse(data2$dement_consensus == "M", NA, data2$status)

#check possible cases were removed
table(data2$dement_consensus) # removed 7 maybes 
#0   1
#452 110

#make list of episcores
list_e <- colnames(data2)[9:92]

#make empty results dataframe
results <- data.frame(SeqId = list_e, Odds_Ratio = NA, N = NA, ci.lower = NA, ci.upper = NA, P = NA, logOdds = NA, SE = NA)

#make episcores rownames
rownames(results) = results$SeqId # This will allow you to index with results[i,]

#list of episcores 
episcores = colnames(data2)[which(colnames(data2)%in% list_e)]

#for loop that loops over episcores and performs logistic regression models
for(i in episcores) { 
  data2$tmp = data2[,i]
  
  mod <- glm(data2$status ~ data2$tmp + scale(data2$age) + as.factor(data2$sex), family= "binomial")
  
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
write.csv(results, file = "LBC1921_dementia_84_episcores_logistic_basic_model_04072023.csv")

####################################################################################################################################################################################################################
####################################################################################################################################################################################################################
#make list of episcores
list_e <- colnames(data2)[9:92]

#make empty results dataframe 
results <- data.frame(SeqId = list_e, Odds_Ratio = NA, N = NA, ci.lower = NA, ci.upper = NA, P = NA, logOdds = NA, SE = NA)

#set episcore names as rownames
rownames(results) = results$SeqId # This will allow you to index with results[i,]

#list of episcores
episcores = colnames(data2)[which(colnames(data2)%in% list_e)]

#for loop that loops over episcores and performs logistic regression models 
for(i in episcores) { 
  data2$tmp = data2[,i]
  
  mod <- glm(data2$status ~ data2$tmp + scale(data2$age) + as.factor(data2$sex) + scale(data2$bmi) + scale(data2$alcpw) + scale(data2$soclcode) + scale(data2$smokingScore), family= "binomial")
  
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

#saves out results
write.csv(results, file = "LBC1921_dementia_84_episcores_logistic_full_model_no_yrsedu_04072023.csv")
