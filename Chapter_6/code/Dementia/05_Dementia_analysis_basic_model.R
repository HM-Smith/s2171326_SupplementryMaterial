library(dplyr)
library(readr)
library(survival)
library(kinship2)
library(coxme)
# 
#read in pedigree file
ped <- read.csv("2023-03-20_pedigree.csv")

#make kinship matrix
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
time_to_event <- read.csv("time_to_event_data_dementia_death_25102024.csv")
#convert to data frame 
data <- as.data.frame(time_to_event)

#create status column where 1 = dementia, 0 = non-affected 
data$status <- ifelse(data$dementia_event == 1, 1, 0)

#create time columm which is either time to dementia, or time to death or time to censor
data$time <- ifelse(data$dementia_event == 1, data$tte_dementia, data$tte_death)



#check for individuasl with negative tte for dementia 
table(data$tte_dementia < 0)
# FALSE  TRUE
# 7127     1


# 1 individual (53681) has dementia diagnosis prior to baseline app

#recode individual to NA 
data$status <- ifelse(data$tte_dementia < 0, NA, data$status)

#check if anyone has negative tte for death 
table(data$tte_death < 0)
# FALSE
# 7128

#check dementia/non-affected N's 
table(data$status)

# 0    1
# 6911  216

#read in protein 
prot <- readRDS("GS_ProteinGroups_RankTransformed_23Aug2023.rds")

names(prot)[1] <- "Sample_Name"

#merge datasets together
data <- left_join(prot, data, by = "Sample_Name")

#read in age and sex data
agesex <- read.csv("2023-03-17_agesex.csv")

names(agesex)[1] <- "Sample_Name"

#join datasets 
data <- left_join(data, agesex, by = "Sample_Name")

#subset to complete data 
data <- data[complete.cases(data), ]

#create list of proteins
list_e <- colnames(data)[c(2:440)]

#create empty results dataframe
results <- data.frame(id = list_e, Outcome = NA, Hazard_Ratio = NA,  ci.lower = NA,ci.upper = NA, P = NA, N_cases = NA, N_controls = NA, local = NA, global = NA)

#make EpiScores rownames
rownames(results) = results$id # This will allow you to index with results[i,]


#for loop for looping over proteins and performing Cox mixed effects models
proteins = colnames(data)[which(colnames(data)%in% list_e)]
for(i in proteins) { 
  data$tmp = data[,i]

 mod = coxme(Surv(time, status) ~  scale(tmp) + factor(sex) + Age_at_baseline + (1|Sample_Name), varlist = kin_model*2, data = data)
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
write.csv(results, file = "Dementia_GS_MS_prots_basic_model_coxme.csv")
