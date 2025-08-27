library(tidyverse)
library(readr)
library(survival)
library(kinship2)
library(coxme)



## Read in pedigree information and make kinship matrix
ped <- read.csv("2023-03-20_pedigree.csv")
kin <- with(ped, pedigree(volid, father, mother, sex, famid=famid))
kin_model <- kinship(kin) 

# Load function to Extract Lmekin Results	
extract_coxme_table <- function (mod){
  beta <- mod$coefficients #$fixed is not needed
  nvar <- length(beta)
  nfrail <- nrow(mod$var) - nvar
  se <- sqrt(diag(mod$var)[nfrail + 1:nvar])
  z<- round(beta/se, 2)
  p<- 1 - pchisq((beta/se)^2, 1)
  table=data.frame(cbind(beta,se,z,p))
  return(table)
}


#load time to event data
time_to_event <- read.csv("time_to_event_data_dementia_death_GenScot_84_EpiScores_15082023.csv")

#convert to data frame 
data <- as.data.frame(time_to_event)

#covert sex to factor 
data$sex <- as.factor(data$sex)

#create status column where 1 = dementia, 0 = non-affected 
data$status <- ifelse(data$dementia_event == 1, 1, 0)

#create time columm which is either time to dementia, or time to death or time to censor
data$time <- ifelse(data$dementia_event == 1, data$tte_dementia, data$tte_death)



#check for individuasl with negative tte for dementia 
table(data$tte_dementia < 0)
#FALSE  TRUE
#7790     1
# 1 individual (53681) has dementia diagnosis prior to baseline app

#recode individual to NA 
data$status <- ifelse(data$tte_dementia < 0, NA, data$status)

#check if anyone has negative tte for death 
table(data$tte_death < 0)
#FALSE
#7791

#check dementia/non-affected N's 
table(data$status)
#0    1
#7555  235

#check dementia/non-affected straified by sex
table(data$status, data$sex) # 0 = males , 1 = females

##  0    1
#0 3096 4459
#1   78  157



#create list of EpiScores
list_e <- colnames(data)[9:92]

#create empty results dataframe
results <- data.frame(SeqId = list_e, Outcome = NA, Hazard_Ratio = NA,  ci.lower = NA,ci.upper = NA, P = NA, N_cases = NA, N_controls = NA, local = NA, global = NA)

#make EpiScores rownames
rownames(results) = results$SeqId # This will allow you to index with results[i,]

#for loop for looping over EpiScores and performing Cox mixed effects models
episcores = colnames(data)[which(colnames(data)%in% list_e)]
for(i in episcores) { 
  data$tmp = data[,i]
  
mod = coxme(Surv(data$time, data$status) ~ data$tmp + factor(data$sex) + scale(data$Age_at_baseline) + (1|data$Sample_Name), varlist = kin_model*2)

results[i,1] <- i
results[i,2] <- as.character("Dementia")
results[i,3:5]<-round(exp(cbind(coef(mod), confint(mod)))[1,1:3],2)
results[i,6] <- extract_coxme_table(mod)[1,4]
results[i,7] <- mod$n[1]
results[i,8] <- mod$n[2]-mod$n[1]

all <- cox.zph(mod) 
p <- all$table[,"p"]
local <- p[1]
global <- p[4]

results[i,9] <- local
results[i,10] <- global

print(i)
}

#save out results dataframe 
write.csv(results, file = "Coxme_Dementia_GenScot_84_episcores_basic_15082023.csv")

###########################################################################################################################################################################################

#read in covariates files
alc <- read.csv("2023-02-01_alc.csv")
edu <- read.csv("2023-02-01_education.csv")
smid <- read.csv("2023-02-13_simd.csv")
BMI <- read.csv("body.csv")

#rename columns to match 
names(alc)[1] <- "Sample_Name"
names(edu)[1] <- "Sample_Name"
names(edu)[2] <- "edu_years"
names(smid)[1] <- "Sample_Name"
names(BMI)[1] <- "Sample_Name"

#select columns of interest
alc <- alc[c(1,2)]
edu <- edu[c(1,2)]
BMI <- BMI[c(1,4)]

#read in and combine EpiSmokEr scores
w1 <- readRDS("wave1_epismoker.rds")
w3 <- readRDS("wave3_epismoker.rds")
w4 <- readRDS("wave4_epismoker.rds")
bind <- rbind(w1, w3)
bind <- rbind(bind, w4) # 18779 individuals with DNAm epismoker calculated 
bind$Sample_Sentrix_ID <- row.names(bind)
bind <- bind[-1]

#join dementia data with covariates
data <- data %>% left_join(alc, by ="Sample_Name")
data <- data %>% left_join(edu, by ="Sample_Name")
data <- data %>% left_join(smid, by ="Sample_Name")
data <- data %>% left_join(BMI, by ="Sample_Name")
data <- data %>% left_join(bind, by ="Sample_Sentrix_ID")

#check for missing values
table(is.na(data))
#FALSE   TRUE
#713478   1522

#removing covariate data points 3.5 SD away from mean  
list <- c("units", "bmi", "rank", "edu_years", "smokingScore")
for(i in list){ 
  
  cutoff1 = mean(data[,i], na.rm = T) + 3.5*sd(data[,i], na.rm = T)
  cutoff2 = mean(data[,i], na.rm = T) - 3.5*sd(data[,i], na.rm = T)
  
  data[,i][which(data[,i] > cutoff1 | data[,i] < cutoff2)] <- NA 
}

#check missing values
table(is.na(data))
#FALSE   TRUE
#713311   1689

#make list of EpiScores
list_e <- colnames(data)[9:92]

#make empty results dataframe
results <- data.frame(SeqId = list_e, Outcome = NA, Hazard_Ratio = NA,  ci.lower = NA,ci.upper = NA, P = NA, N_cases = NA, N_controls = NA, local = NA, global = NA)

#set EpiScores as rownames
rownames(results) = results$SeqId # This will allow you to index with results[i,]

#loop for looping over EpiScores and performing Cox mixed effects models 
episcores = colnames(data)[which(colnames(data)%in% list_e)]
for(i in episcores) { 
  data$tmp = data[,i]
  
  mod = coxme(Surv(data$time, data$status) ~ scale(data$tmp) + factor(data$sex) + scale(data$Age_at_baseline) + scale(data$units) + scale(data$bmi) + scale(data$rank)  + scale(data$smokingScore) + (1|data$Sample_Name), varlist = kin_model*2)
  
  results[i,1] <- i
  results[i,2] <- as.character("Dementia")
  results[i,3:5]<-round(exp(cbind(coef(mod), confint(mod)))[1,1:3],2)
  results[i,6] <- extract_coxme_table(mod)[1,4]
  results[i,7] <- mod$n[1]
  results[i,8] <- mod$n[2]-mod$n[1]
  
  all <- cox.zph(mod) 
  p <- all$table[,"p"]
  local <- p[1]
  global <- p[8]
  
  results[i,9] <- local
  results[i,10] <- global
  
  print(i)
}

#save out results
write.csv(results, file = "Coxme_Dementia_GenScot_84_episcores_Full_model_no_yrsedu_15082023.csv")


