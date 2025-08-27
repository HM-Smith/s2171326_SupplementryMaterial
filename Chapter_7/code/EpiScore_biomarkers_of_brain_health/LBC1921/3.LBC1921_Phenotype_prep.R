## Data in Spss file, converted to csv file
library("haven")
df <- read_sav("LBC1921_HannahSmith_EpigeneticProteinScoresAndBrainHealth_24JAN2023.sav")
df <- as.data.frame(df)
write.csv(df, "LBC1921_HannahSmith_EpigeneticProteinScoresAndBrainHealth_24JAN2023.csv")

#load phenotype data 
library(readr)
LBC1921 <- read.csv("LBC1921_HannahSmith_EpigeneticProteinScoresAndBrainHealth_24JAN2023.csv")

## Data inspection 
dim(LBC1921) 
##[1] 569 102


#empty lists
na_neg888 = na_999 = na_neg999 = na_neg99 = na_neg777 = list()

#check NAs
table(is.na(LBC1921))
##FALSE  TRUE
##26117 31921


#loop over columns and convert missing codes to NA, and record how many were changed
for(i in 4:ncol(LBC1921)){
  na_neg888[[i]] = length(which(as.numeric(LBC1921[,i]) == -888))
  na_999[[i]] = length(which(as.numeric(LBC1921[,i]) == 999))
  na_neg999[[i]] = length(which(as.numeric(LBC1921[,i]) == -999))
  na_neg99[[i]] = length(which(as.numeric(LBC1921[,i]) == -99))
  na_neg777[[i]] = length(which(as.numeric(LBC1921[,i]) == -777))
  
  LBC1921[,i][which(as.numeric(LBC1921[,i]) == -999)] = NA 
  LBC1921[,i][which(as.numeric(LBC1921[,i]) == -888)] = NA 
  LBC1921[,i][which(as.numeric(LBC1921[,i]) == -777)] = NA 
  LBC1921[,i][which(as.numeric(LBC1921[,i]) == -99)] = NA 
  LBC1921[,i][which(as.numeric(LBC1921[,i]) == 999)] = NA 
}

# check how many NA's after recoding
table(is.na(LBC1921))
#FALSE  TRUE
#26107 31931

## how many NA's introduced by recoding of SPSS codes 
31931 - 31921 # [1] 10


## does this match list of SPSS codes
a = unlist(na_neg777) + unlist(na_neg888) + unlist(na_neg99) + unlist(na_neg999) + unlist(na_999)
sum(a) #[1] 10


## how many of each code
library(tidyverse)
unlist(na_neg777) %>% sum #0
unlist(na_neg888) %>% sum #0 
unlist(na_neg99) %>% sum #10
unlist(na_neg999) %>% sum #0
unlist(na_999) %>% sum #0



## recode sex to 0 and 1
library(tidyverse)
LBC1921 <- mutate(LBC1921, sex=recode(gender, '1' = '0', '2' = '1')) # 0 = male, 1 = female 
#convert sex to factor variable 
LBC1921$sex <- as.factor(LBC1921$sex)


# check dataset
summary(LBC1921)

## read in LBC1921 episcores 
LBC1921_EpiScores <- read.csv("KORA_EpiScores_projected_Adjusted_in_LBC1921_WAVE1_29012023.csv")

## rename ID to lbc36no
LBC1921_EpiScores <- LBC1921_EpiScores %>% rename(studyno = ID)

## select columns 
LBC1921_EpiScores = LBC1921_EpiScores %>% select(3:ncol(LBC1921_EpiScores))
## left join to LBC1921 dataset 
LBC1921_epi <- LBC1921 %>% left_join(LBC1921_EpiScores, by = "studyno")

#rename sex.x and set as factor 
LBC1921_epi <- LBC1921_epi %>% rename(sex = sex.x)
LBC1921_epi$sex <- as.factor(LBC1921_epi$sex)

# write out file 
write.csv(LBC1921_epi, file = "LBC1921_WAVE1_cohort_data_KORA_EpiScore_Adjusted_31012023.csv")
