library(metafor)
library(dplyr)

#read in dementia diagnosis result files 
genscot_dem <- read.csv("logistic_Dementia_GenScot_84_episcores_with_ids_basic_16082023.csv")

LBC1921_dem <- read.csv("LBC1921_dementia_84_episcores_with_ids_logistic_basic_model_05072023.csv")

LBC1936_dem <- read.csv("logistic_Dementia_LBC1936_84_EpiScores_basic_model_WITH_IDS_20072023.csv")

#rename column to match
names(genscot_dem)[4]<- "Odds_Ratio"

#select columns of interest 
genscot_dem <- genscot_dem %>% select(3:14)
LBC1921_dem <- LBC1921_dem %>% select(3:14)
LBC1936_dem <- LBC1936_dem %>% select(3:14)

#make list of episcores
EpiScores <- genscot_dem$SeqId 

#set episcores as rownames of results files 
rownames(genscot_dem) <- genscot_dem$SeqId
rownames(LBC1921_dem) <- LBC1921_dem$SeqId
rownames(LBC1936_dem) <- LBC1936_dem$SeqId

#create empty results dataframe
outcomes <- data.frame("SeqId" = EpiScores, "Gene_Name" = NA, "Target" = NA, "Uniprot"= NA, "q" = NA, "p_for_q" = NA, "I^2" = NA, "tau^2" = NA, "setau^2" = NA, "beta" = NA, 
                       "se" = NA, "z" = NA, "p" = NA, "H^2" = NA, "lower_CI" = NA, "Upper_CI" = NA)

#set episcores as rownames
rownames(outcomes) <- outcomes$SeqId

#for loop for looping over episcores and performing meta-analysis 
for (i in EpiScores) {
  loopdata <- rbind(genscot_dem[i,], LBC1921_dem[i,], LBC1936_dem[i,]) # dataframe with the episcores per row and columns with betas and variances, and can loop over the rows
  
  
  res <- rma(logOdds, sei = SE, data = loopdata) # meta-analysis
  
  outcomes[i,1] <- loopdata$SeqId[1]
  outcomes[i,2] <- loopdata$Gene_Name[1]
  outcomes[i,3] <- loopdata$Target[1]
  outcomes[i,4] <- loopdata$UniProt[1]
  
  outcomes[i,5] <- as.numeric(res[17]) # q statistic 
  outcomes[i,6] <- as.numeric(res[18]) # p value for q
  outcomes[i,7] <- as.numeric(res[13]) # I^2
  outcomes[i,8] <- as.numeric(res[9]) # tau^2
  outcomes[i,9] <- as.numeric(res[10]) # setau^2
  outcomes[i,10] <- as.numeric(res[2]) #beta
  outcomes[i,11] <- as.numeric(res[3]) # se
  outcomes[i,12] <- as.numeric(res[4]) # z
  outcomes[i,13] <- as.numeric(res[5]) # p value
  outcomes[i,14] <- as.numeric(res[14]) # H^2
  outcomes[i,15] <- as.numeric(res[6])# lower CI
  outcomes[i,16] <- as.numeric(res[7])# upper CI 
  
}

#create FDR corrected column 
outcomes$FDR <- p.adjust(outcomes$p, method = "BH")

#save out all results
write.csv(outcomes, file = "dementia_logistic(binary)_all_three_cohorts_meta_analysis_basic_results_16082023.csv")

#filter to FDR signficant results
outcomes_sig <- outcomes %>% filter(FDR < 0.05)

#check number of FDR significant results
dim(outcomes_sig)
#[1] 0 17

########################################################################################################################################################################################################################################################
########################################################################################################################################################################################################################################################
library(metafor)
library(dplyr)

genscot_dem <- read.csv("logistic_Dementia_GenScot_84_episcores_with_ids_Full_model_no_yrsedu_16082023.csv")

LBC1921_dem <- read.csv("LBC1921_dementia_84_episcores_with_ids_logistic_full_model_no_yrsedu_05072023.csv")

LBC1936_dem <- read.csv("logistic_Dementia_LBC1936_84_EpiScores_full_model_no_yrsedu_WITH_IDS_27102023.csv")

#select columns of interest 
genscot_dem <- genscot_dem %>% select(3:14)
LBC1921_dem <- LBC1921_dem %>% select(3:14)
LBC1936_dem <- LBC1936_dem %>% select(3:14)

#make list of episcores
EpiScores <- genscot_dem$SeqId 

#set episcores as rownames of results files
rownames(genscot_dem) <- genscot_dem$SeqId
rownames(LBC1921_dem) <- LBC1921_dem$SeqId
rownames(LBC1936_dem) <- LBC1936_dem$SeqId

#make empty results dataframe
outcomes <- data.frame("SeqId" = EpiScores, "Gene_Name" = NA, "Target" = NA, "Uniprot"= NA, "q" = NA, "p_for_q" = NA, "I^2" = NA, "tau^2" = NA, "setau^2" = NA, "beta" = NA, 
                       "se" = NA, "z" = NA, "p" = NA, "H^2" = NA, "lower_CI" = NA, "Upper_CI" = NA)

#set episcores as rownames 
rownames(outcomes) <- outcomes$SeqId


#for loop for looping over episcores and performing meta-analysis 
for (i in EpiScores) {
  loopdata <- rbind(genscot_dem[i,], LBC1921_dem[i,], LBC1936_dem[i,]) # dataframe with the episcores per row and columns with betas and variances, and can loop over the rows
  
  
  res <- rma(logOdds, sei = SE, data = loopdata) # meta-analysis
  
  outcomes[i,1] <- loopdata$SeqId[1]
  outcomes[i,2] <- loopdata$Gene_Name[1]
  outcomes[i,3] <- loopdata$Target[1]
  outcomes[i,4] <- loopdata$UniProt[1]
  
  outcomes[i,5] <- as.numeric(res[17]) # q statistic 
  outcomes[i,6] <- as.numeric(res[18]) # p value for q
  outcomes[i,7] <- as.numeric(res[13]) # I^2
  outcomes[i,8] <- as.numeric(res[9]) # tau^2
  outcomes[i,9] <- as.numeric(res[10]) # setau^2
  outcomes[i,10] <- as.numeric(res[2]) #beta
  outcomes[i,11] <- as.numeric(res[3]) # se
  outcomes[i,12] <- as.numeric(res[4]) # z
  outcomes[i,13] <- as.numeric(res[5]) # p value
  outcomes[i,14] <- as.numeric(res[14]) # H^2
  outcomes[i,15] <- as.numeric(res[6])# lower CI
  outcomes[i,16] <- as.numeric(res[7])# upper CI 
  
}

#create FDR corrected column 
outcomes$FDR <- p.adjust(outcomes$p, method = "BH")

#save out all results 
write.csv(outcomes, file = "dementia_logistic(binary)_all_three_cohorts_meta_analysis_full_model_no_years_edu_results_30102023.csv")

#filter to FDR significant results
outcomes_sig <- outcomes %>% filter(FDR < 0.05)

#check number of FDR significant results 
dim(outcomes_sig)
#[1] 0 17

