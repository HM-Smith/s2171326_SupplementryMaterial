####### LBC1936 PROJECT ###########
## Data in Spss file, converted to csv file
library("haven")
df <- read_sav("LBC1936_HannahSmithEpigeneticProteinScoresAndBrainHealth_SH_04OCT2022.sav")
df <- as.data.frame(df)
write.csv(df, "LBC1936_HannahSmith_19122022.csv")

#load phenotype data
library(tidyverse)
LBC1936 <- read.csv("LBC1936_HannahSmith_19122022.csv")

## Data inspection 
dim(LBC1936) #1091 rows & 326 columns

#make lists to track missing codes 
na_neg888 = na_999 = na_neg999 = na_neg99 = na_neg777 = list()

#check NA's 
table(is.na(LBC1936))
# FALSE   TRUE 
# 209137 146529 


#loop that loops over columns and recodes missing codes to NA and tracks how many of each code
for(i in 3:ncol(LBC1936)){
  na_neg888[[i]] = length(which(as.numeric(LBC1936[,i]) == -888))
  na_999[[i]] = length(which(as.numeric(LBC1936[,i]) == 999))
  na_neg999[[i]] = length(which(as.numeric(LBC1936[,i]) == -999))
  na_neg99[[i]] = length(which(as.numeric(LBC1936[,i]) == -99))
  na_neg777[[i]] = length(which(as.numeric(LBC1936[,i]) == -777))
  
  LBC1936[,i][which(as.numeric(LBC1936[,i]) == -999)] = NA 
  LBC1936[,i][which(as.numeric(LBC1936[,i]) == -888)] = NA 
  LBC1936[,i][which(as.numeric(LBC1936[,i]) == -777)] = NA 
  LBC1936[,i][which(as.numeric(LBC1936[,i]) == -99)] = NA 
  LBC1936[,i][which(as.numeric(LBC1936[,i]) == 999)] = NA 
}

# check how many NA's after recoding
table(is.na(LBC1936))
#FALSE   TRUE 
#207279 148387

## how many NA's introduced by recoding of SPSS codes 
148387 - 146529 # [1] 1858

## does this match list of SPSS codes
a = unlist(na_neg777) + unlist(na_neg888) + unlist(na_neg99) + unlist(na_neg999) + unlist(na_999)
sum(a) # 1858

## how many of each code
unlist(na_neg777) %>% sum #4
unlist(na_neg888) %>% sum # 16
unlist(na_neg99) %>% sum #296
unlist(na_neg999) %>% sum #1531
unlist(na_999) %>% sum #11



## recode sex to 0 and 1
LBC1936 <- mutate(LBC1936, sex=recode(sex, '1' = '0', '2' = '1')) # 0 = male, 1 = female 
#check sex to factor variable 
LBC1936$sex <- as.factor(LBC1936$sex)

## age in days converted to years
LBC1936$Age_w1 <- LBC1936$agedays_w1 / 365.25
LBC1936$Age_w2 <- LBC1936$agedays_w2 / 365.25
LBC1936$Age_w3 <- LBC1936$agedays_w3 / 365.25
LBC1936$Age_w4 <- LBC1936$agedays_w4 / 365.25
LBC1936$Age_w5 <- LBC1936$agedays_w5 / 365.25


# check dataset
summary(LBC1936)

## recoding impossible values to NA
check_symsear_w1 <- LBC1936 %>% filter(symsear_w1 < 0) %>% select(symsear_w1) # one person has -1, one person has -4 
check_symsear_w4 <- LBC1936 %>% filter(symsear_w4 <0) %>% select(symsear_w4) # one person has -3 

LBC1936$symsear_w1 =ifelse(LBC1936$symsear_w1 < 0, NA, LBC1936$symsear_w1)
LBC1936$symsear_w4 =ifelse(LBC1936$symsear_w4 %in% -3, NA, LBC1936$symsear_w4)


## read in LBC1936 episcores 
LBC1936_EpiScores <- read.csv("KORA_EpiScores_projected_Adjusted_in_LBC1936_WAVE1_17072023.csv")

## rename ID to lbc36no
LBC1936_EpiScores <- LBC1936_EpiScores %>% rename(lbc36no = ID)

## left join to LBC1936 dataset 
LBC1936_epi <- LBC1936 %>% left_join(LBC1936_EpiScores, by = "lbc36no")

#rename sex.x and set as factor 
LBC1936_epi <- LBC1936_epi %>% rename(sex = sex.x)
LBC1936_epi$sex <- as.factor(LBC1936_epi$sex)

# write out file 
write.csv(LBC1936_epi, file = "LBC1936_WAVE1_cohort_data_KORA_EpiScore_Adjusted_17072023.csv")
