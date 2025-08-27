# Load in packages 
library(survival)
library(coxme)
library(readxl)
library(tidyverse)
library(gdata)
library(haven)
library(foreign)


# Load dementia diagnoses
dem <- read_sav("LBC1936_EarlyAccessDementiaAscertainment.sav")

#load in death info
death <- read_sav("LBC1936_HannahSmithEpigeneticProteinScoresAndBrainHealth_AgeDeath_17FEB2023.sav")

#load in screening info 
screening <- read_sav("DemAcertainmentScreeningAge.sav")

#set as dataframes
death <- as.data.frame(death)
dem <- as.data.frame(dem)
screening <- as.data.frame(screening)


# How many cases of dementia
table(dem$ConsensusDementiaDiagnosis)

#              -888              -999          Dementia               MCI
#               702               226               118                13
#       No dementia       No Dementia      Not dementia Possible dementia
#                 7                 2                16                 7

table(dem$dementia_code)
#   0   1 
# 747 118



# Filter to W2 LBC1936 target and select columns of interest
target <- readRDS("targets_3489_bloodonly.rds")
dat3 <- target[which(target$WAVE %in% '2'),]
dat3  <- dat3 %>% filter(cohort == "LBC36") 
dim(dat3)# 792 people 
dat3 <- dat3 %>% select('Basename', 'set', 'ID')


## check dementia file 
dim(dem)
#[1] 1091    6

# Subset dementia file to W2 individuals
sub <- dem[which(dem$lbc36no %in% dat3$ID),]
#

#check subset data 
dim(sub)
#[1] 792   6


# Look at dementia diagnoses in the context of W2

table(sub$ConsensusDementiaDiagnosis)

#-888              -999          Dementia               MCI
#641                 1               107                13
#No dementia       No Dementia      Not dementia Possible dementia
#6                 2                15                 7

   
table(sub$dementia_code)

#   0   1
# 684 107

# Remove individuals with possible dementia or MCI or missing codes
sub <- sub[-which(sub$ConsensusDementiaDiagnosis %in% c('MCI', 'Possible dementia', '-999')),] # 780

#check subset file
dim(sub)
#[1] 780   6

#check dementia cases 
table(sub$ConsensusDementiaDiagnosis)
# -888     Dementia  No dementia  No Dementia Not dementia
#  649          108            6            2           15

#check dementia cases 
table(sub$dementia_code)
#   0   1
# 671 108


# 108 diagnoses and assign all those without dementia as controls
# Integrate DOB imputed to July 1936, with censor date set as per screening file 
# Age event will be calculated depending on case/control and alive/dead status in the study window

#rename id column 
names(dat3)[3] <- 'lbc36no'

####################################################################################

### INCDIENT ANALYSES

####################################################################################

#### PREP PHENOTYPE FILE WITH AGE ALIVE AND AGE DEATH INFO 

# Calculate age at censor date for all individuals (dead or alive)
d1 <- sub %>% left_join(screening, by = "lbc36no")
d1$ageyears_censor <- (d1$agedays_DemAsc_w6 / 365.25)

#check d1 file
dim(d1) 
#[1] 771  8

#combined d1 and death info
d1 <- d1 %>% left_join(death, by = "lbc36no")

#check d1 file 
dim(d1)
#[1] 771  10

#remove NA's from age at death 
dead <- d1[complete.cases(d1$ageyrs_death),] 

#check number of deaths 
dim(dead)
#[1] 285  10


# Map those that died before the censor date (i.e. age dead occurs before age at censor)
dead$died_pre_censor <- ifelse(dead$ageyrs_death < dead$ageyears_censor, 1, 0) 
table(dead$died_pre_censor) # all died pre censor date 
#0   1
#8 277


# Assign dead status to main file
d1$dead <- ifelse(d1$ageyrs_death > 0, 1, 0)
d1$dead[is.na(d1$dead)] <- 0
table(d1$dead)
#0   1
#486 2


# Create variable with age dead if died, or age at censor if alive at censor 
d1$aged <- ifelse(d1$dead == '1' & d1$ageyrs_death < d1$ageyears_censor, d1$ageyrs_death, d1$ageyears_censor)

# read in LBC1936 data
LBC1936 <- read.csv("LBC1936_WAVE1_cohort_data_KORA_EpiScore_Adjusted_17072023.csv")

#select columns of interest
LBC1936_age <- LBC1936 %>% select(lbc36no, agedays_w1, sex, alcunitwk_w1, bmi_w1, yrsedu_w1, depind_w1)

#join d1 and LBC36 data
d1 <- d1 %>% left_join(LBC1936_age, by = "lbc36no")

#convert age days to years
d1$age_baseline = d1$agedays_w1/365.25

#calculate time to death from baseline 
d1$tte_death <- as.numeric(d1$aged - d1$age_baseline)

#set dementia code as numeric
d1$dementia_code <- as.numeric(d1$dementia_code)

#create age at dementia column 
d1$age_at_dementia <- ifelse(d1$dementia_code == "1", d1$age_at_dementia_diagnosis, d1$aged)

#calculate time to dementia column 
d1$tte_dementia <- as.numeric(d1$age_at_dementia - d1$age_baseline)

#check if any negative time to dementia 
table(d1$tte_dementia < 0)

#FALSE
#771

 #check if any negative time to death   
table(d1$tte_death < 0)

#FALSE
#771

#check if any time to dementia is after time of death 
table(d1$tte_dementia > d1$tte_death)

#FALSE  TRUE
#768     3

#check if time to dementia equals time to death 
table(d1$tte_dementia == d1$tte_death)

#FALSE  TRUE
#107   664

#check dementia cases by death cases 
table(d1$dementia_code, d1$dead)

#0   1
#0 438 226
#1  48  59

#set dementia and death code as numeric 
d1$dementia_code <- as.numeric(d1$dementia_code)
d1$dead <- as.numeric(d1$dead)


summary(d1$age_at_dementia_diagnosis) # 796 included still 

#check d1 file 
dim(d1)
#[1] 771  22

   

# Add episcore markers that have been pre-corrected for DNAm technical variables and rank-inverse based transformed
scores <- read.csv("LBC1936_WAVE1_cohort_data_KORA_EpiScore_Adjusted_17072023.csv")
scores <- scores %>% select(lbc36no,335:419)

# join episcores with dementia/death data 
cox <- left_join(d1, scores, by = 'lbc36no')

#read in and join epismoker
epismoke <- readRDS("lbc_epismoker.rds") 
names(epismoke)[1] <- "Basename"
cox <- cox %>% left_join(epismoke, by = "Basename")

dim(cox)
#[1] 771 108

   
#save out file 
write.csv(cox, file = "time_to_event_dementia_death_data_LBC1936_new_ascertainment_20072023.csv")

