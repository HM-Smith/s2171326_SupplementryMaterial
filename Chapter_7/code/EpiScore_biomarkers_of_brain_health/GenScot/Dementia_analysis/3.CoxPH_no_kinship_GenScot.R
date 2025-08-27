library(tidyverse)
library(readr)
library(survival)
library(kinship2)
library(coxme)

#read in time to event data
data <- read.csv("time_to_event_data_dementia_death_GenScot_84_EpiScores_15082023.csv")

#set sex as a factor
data$sex <- as.factor(data$sex)

#create status column where dementia = 1 and un-affected = 0
data$status <- ifelse(data$dementia_event == 1, 1, 0)

#create time column with either time to dementia or time to death or time to censor 
data$time <- ifelse(data$dementia_event == 1, data$tte_dementia, data$tte_death)

#remove anyone diagnosed before baseline 
data$status <- ifelse(data$tte_dementia < 0, NA, data$status)

#create list of EpiScores
list_e <- colnames(data)[9:92]

#create empty results dataframe
results <- data.frame(SeqId = list_e, Outcome = NA, Hazard_Ratio = NA,  ci.lower = NA,ci.upper = NA, P = NA, N_cases = NA, N_controls = NA, local = NA, global = NA)

#set episcores are rownames
rownames(results) = results$SeqId # This will allow you to index with results[i,]

#for loop that loops over episcores and performs CoxPH models 
episcores = colnames(data)[which(colnames(data)%in% list_e)]
for(i in episcores) { 
  data$tmp = data[,i]
  
  mod = coxph(Surv(data$time, data$status) ~ data$tmp + factor(data$sex) + data$Age_at_baseline)
  
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

#save out results file
write.csv(results, file = "coxph_no_kinship_GenScot_dementia_basic_model_15082023.csv", row.names = F)

###########################################################################################################################################################################################
#read in covariates files
alc <- read.csv("2023-02-01_alc.csv")
smid <- read.csv("2023-02-13_simd.csv")
BMI <- read.csv("body.csv")

#rename columns 
names(alc)[1] <- "Sample_Name"
names(smid)[1] <- "Sample_Name"
names(BMI)[1] <- "Sample_Name"

#select columns
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

#combine dementia data and covariates 
data <- data %>% left_join(alc, by ="Sample_Name")
data <- data %>% left_join(smid, by ="Sample_Name")
data <- data %>% left_join(BMI, by ="Sample_Name")
data <- data %>% left_join(bind, by ="Sample_Sentrix_ID")

#check missing values
table(is.na(data))
#FALSE   TRUE
#713478   1522

#remove covariate data points 3.5 SD from mean 
list <- c("units", "bmi", "rank", "smokingScore")
for(i in list){ 
  
  cutoff1 = mean(data[,i], na.rm = T) + 3.5*sd(data[,i], na.rm = T)
  cutoff2 = mean(data[,i], na.rm = T) - 3.5*sd(data[,i], na.rm = T)
  
  data[,i][which(data[,i] > cutoff1 | data[,i] < cutoff2)] <- NA 
}

#check missing values
table(is.na(data))
#FALSE   TRUE
#713311   1689

#make list of episcores
list_e <- colnames(data)[9:92]

#make empty resuts dataframe
results <- data.frame(SeqId = list_e, Outcome = NA, Hazard_Ratio = NA,  ci.lower = NA,ci.upper = NA, P = NA, N_cases = NA, N_controls = NA, local = NA, global = NA)

#set rownames as episcores
rownames(results) = results$SeqId # This will allow you to index with results[i,]

#for loop that loops over episcores and performs coxPH models 
episcores = colnames(data)[which(colnames(data)%in% list_e)]
for(i in episcores) { 
  data$tmp = data[,i]
  
  mod = coxph(Surv(data$time, data$status) ~ scale(data$tmp) + factor(data$sex) + scale(data$Age_at_baseline) + scale(data$units) + scale(data$bmi) + scale(data$rank)  + scale(data$smokingScore))
  
  mod.sum <- summary(mod)$coefficients
  mod.conf <- confint(mod)
  results[i,1] <- i
  results[i,2] <- as.character("Dementia")
  results[i,3]<- mod.sum[1,2]
  results[i,4]<- exp(mod.conf[1,1])
  results[i,5]<- exp(mod.conf[1,2])
  
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
write.csv(results, file = "CoxPH_GenScot_dementia_84_episcores_Full_model_no_yrsedu_15082023.csv")


