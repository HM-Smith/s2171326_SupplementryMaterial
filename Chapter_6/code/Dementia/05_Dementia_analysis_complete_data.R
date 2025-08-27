library(dplyr)
library(readr)
library(survival)
library(kinship2)
library(coxme)
# 
#read in pedigree file
ped1 <- read.csv("2023-03-20_pedigree.csv")
#input missing fam IDs
ped2 = data.frame(famid=c(4091,4384), volid=c(103027, 144865), father=c(0,0), mother=c(0,0), sex= c("F", "F"))

#bind together
ped = rbind(ped1, ped2)
# 
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


data <- left_join(prot, data, by = "Sample_Name")


covars <- read.csv("covars_lancet_post_outlier_removal.csv")

names(covars)[1] <- "Sample_Name"

data <- left_join(data, covars, by = "Sample_Name")

apoe <- read.csv("apoe_status.csv")

data <- left_join(data, apoe, by = "Sample_Name")

data <- data[complete.cases(data), ]

dim(data)
#[1] 5388  460


table(data$status)
# 0    1
# 5247  141


# Create list of proteins
list_e <- colnames(data)[c(2:440)]

# Create empty results dataframe
results <- data.frame(
  id = list_e,
  Outcome = NA,
  Hazard_Ratio = NA,
  ci.lower = NA,
  ci.upper = NA,
  P = NA,
  N_cases = NA,
  N_controls = NA,
  local = NA,
  global = NA
)

# Make EpiScores rownames
rownames(results) <- results$id  # This will allow you to index with results[i,]

# Loop over proteins and perform Cox mixed effects models
proteins <- colnames(data)[which(colnames(data) %in% list_e)]

# Empty list for storing error proteins
error_proteins <- c()


#for loop for looping over proteins and performing Cox mixed effects models
for (i in proteins) {
  tryCatch({
    data$tmp <- data[, i]
    
    mod = coxme(Surv(time, status) ~  scale(tmp) + factor(sex) + Age_at_baseline + scale(rank) + scale(bmi) + scale(units) + factor(depression_Y) +
                  scale(years) + factor(high_BP_Y) + factor(e4_count) + 
                  (1|Sample_Name), varlist = kin_model*2, data = data)
    results[i,1] <- i
    results[i,2] <- as.character("Dementia")
    results[i,3:5]<-round(exp(cbind(coef(mod), confint(mod)))[1,1:3],2)
    results[i,6] <- extract_coxme_table(mod)[1,4]
    results[i,7] <- mod$n[1]
    results[i,8] <- mod$n[2]-mod$n[1]
    
    all <- cox.zph(mod)
    p <- all$table[,"p"]
    local <- p[1]
    global <- p[11]
    
    results[i,9] <- local
    results[i,10] <- global
    
    # Print progress
    print(i)
  }, error = function(e) {
    # If an error occurs, store error details and skip to next protein
    message(sprintf("Error in processing protein '%s': %s", i, e$message))
    error_proteins <<- c(error_proteins, i)  # Add protein to error list
  })
}
#save out results dataframe 
write.csv(results, file = "Dementia_GS_MS_prots_complete_data_model.csv")
