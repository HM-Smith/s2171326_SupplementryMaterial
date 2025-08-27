library(tidyverse)

#load dementia data 
data <- read.csv("time_to_event_data_dementia_death_GenScot_84_EpiScores_15082023.csv")

#check number of dementia cases, dementia = 1
table(data$dementia_event)

#  0    1
#  7555  236

#select columns of interest 
data2  <- data %>% select(2:3, 7:93)

#read in covariate files
alc <- read.csv("2023-02-01_alc.csv")
smid <- read.csv("2023-02-13_simd.csv")
BMI <- read.csv("body.csv")

#rename columns to match
names(alc)[1] <- "Sample_Name"
names(smid)[1] <- "Sample_Name"
names(BMI)[1] <- "Sample_Name"

#select columns of interest
alc <- alc[c(1,2)]
BMI <- BMI[c(1,4)]

#read in EpiSmokEr scores and combine
w1 <- readRDS("wave1_epismoker.rds")
w3 <- readRDS("wave3_epismoker.rds")
w4 <- readRDS("wave4_epismoker.rds")
bind <- rbind(w1, w3)
bind <- rbind(bind, w4) # 18779 individuals with DNAm epismoker calculated 
bind$Sample_Sentrix_ID <- row.names(bind)
bind <- bind[-1]

#combine dementia and covariates
data2 <- data2 %>% left_join(alc, by ="Sample_Name")
data2 <- data2 %>% left_join(smid, by ="Sample_Name")
data2 <- data2 %>% left_join(BMI, by ="Sample_Name")
data2 <- data2 %>% left_join(bind, by ="Sample_Sentrix_ID")

#check missing values
table(is.na(data2))
#FALSE   TRUE
#670579   1521

#remove covariate data points 3.5 SD from the mean 
list <- c("units", "bmi", "rank", "smokingScore")
for(i in list){ 
  
  cutoff1 = mean(data2[,i], na.rm = T) + 3.5*sd(data2[,i], na.rm = T)
  cutoff2 = mean(data2[,i], na.rm = T) - 3.5*sd(data2[,i], na.rm = T)
  
  data2[,i][which(data2[,i] > cutoff1 | data2[,i] < cutoff2)] <- NA 
}

#check missing values
table(is.na(data2))
#FALSE   TRUE
#670412   1688

#make episcore list
list_e <- colnames(data2)[5:88]

#make empty results dataframe
results <- data.frame(SeqId = list_e, Hazard_Ratio = NA, N = NA, ci.lower = NA, ci.upper = NA, P = NA, logOdds = NA, SE = NA)

#set episcores as rownames
rownames(results) = results$SeqId # This will allow you to index with results[i,]

#for loop that loops over episcores and performs logisitc regression models 
episcores = colnames(data2)[which(colnames(data2)%in% list_e)]
for(i in episcores) { 
  data2$tmp = data2[,i]
  
  mod <- glm(data2$dementia_event ~ data2$tmp + scale(data2$Age_at_baseline) + as.factor(data2$sex), family= "binomial")
  
  Beta <- exp(mod$coefficients['data2$tmp'])
  P <- summary(mod)$coefficients['data2$tmp','Pr(>|z|)']
  CI.lower <- exp(confint(mod)['data2$tmp', 1])
  CI.upper <- exp(confint(mod)['data2$tmp', 2])
  N <- nobs(mod)
  log_HR <- summary(mod)$coefficients["data2$tmp", 1]
  SE <- summary(mod)$coefficients["data2$tmp", 2]
  
  results[i,1] <- i
  results[i,2] <- Beta
  results[i,3] <- N
  results[i,4] <- CI.lower
  results[i,5] <- CI.upper
  results[i,6] <- P
  results[i,7] <- log_HR
  results[i,8] <- SE
  
  print(i)
}

#save out results
write.csv(results, "Logistic_Dementia_GenScot_assocs_84_EpiScores_basic_model_updated_15082023.csv")

####################################################################################################################################################################################################################
####################################################################################################################################################################################################################

#make list of episcores
list_e <- colnames(data2)[5:88]

#create empty results dataframe 
results <- data.frame(SeqId = list_e, Odds_Ratio = NA, N = NA, ci.lower = NA, ci.upper = NA, P = NA, logOdds = NA, SE = NA)

#set episcores as rownames
rownames(results) = results$SeqId # This will allow you to index with results[i,]

#for loop that loops over episcores and performs logistic regression models 
episcores = colnames(data2)[which(colnames(data2)%in% list_e)]
for(i in episcores) { 
  data2$tmp = data2[,i]
  
  mod <- glm(data2$dementia_event ~ data2$tmp + scale(data2$Age_at_baseline) + as.factor(data2$sex) + scale(data2$bmi) + scale(data2$units) + scale(data2$rank) + scale(data2$smokingScore), family= "binomial")
  
  Beta <- exp(mod$coefficients['data2$tmp'])
  P <- summary(mod)$coefficients['data2$tmp','Pr(>|z|)']
  CI.lower <- exp(confint(mod)['data2$tmp', 1])
  CI.upper <- exp(confint(mod)['data2$tmp', 2])
  N <- nobs(mod)
  log_HR <- summary(mod)$coefficients["data2$tmp", 1]
  SE <- summary(mod)$coefficients["data2$tmp", 2]
  
  
  results[i,1] <- i
  results[i,2] <- Beta
  results[i,3] <- N
  results[i,4] <- CI.lower
  results[i,5] <- CI.upper
  results[i,6] <- P
  results[i,7] <- log_HR
  results[i,8] <- SE
  
  print(i)
}

#save out reesults file 
write.csv(results, "Logistic_Dementia_GenScot_assocs_84_EpiScores_full_model_no_yrsedu_updated_15082023.csv")
