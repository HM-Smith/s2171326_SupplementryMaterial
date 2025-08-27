library(tidyverse)


############################################################ CALCULATING Age at Baseline and Age at Oct2020 ####################################################################################################################################

#read in baseline appointment date
baseline_app <- read.csv("2023-02-03_baseline_appt.csv")

#split date into month and year 
baseline_app$baseline_app_year<- as.numeric(substr(baseline_app$ym, 1, 4))
baseline_app$baseline_app_month <- as.numeric(substr(baseline_app$ym, 5, 6))

#read in date of birth file 
dob <- read.csv("2023-02-03_age_yobmob.csv")

#check how many missing month of birth
table(is.na(dob$mob))
#FALSE  TRUE
#22712  1373

table(is.na(dob$yob))
#FALSE
#24085


# Now impute missing MOB information to July (07) as midpoint in year 
dob$mob <- ifelse(is.na(dob$mob), 7, dob$mob)

#check the months have been imputed
table(is.na(dob$mob))
#FALSE
#24085

#combine basline appointment file and dob file 
baseline_app <- baseline_app %>% left_join(dob, by = "id")

#check for missing yob info after merge
table(is.na(baseline_app$yob))
#FALSE  TRUE
#24080    11

#check for missing mob info after merge
table(is.na(baseline_app$mob))
#FALSE  TRUE
#24080    11

###After checking a further age and sex file provided it was confirmed the 
### 11 individuals do not have age,sex or dob info 

 
#calculate age (year) at baseline
baseline_app$year_baseline <- baseline_app$baseline_app_year - baseline_app$yob
#calculate age (month) at baseline 
baseline_app$month_baseline <- (baseline_app$baseline_app_month - baseline_app$mob)/12
#add year and month together to get age age at baseline in years and months
baseline_app$Age_at_baseline <- baseline_app$year_baseline + baseline_app$month_baseline

#calculate age (year) at censor (Apr 2022)
baseline_app$diff_to_2022 <- 2022 - baseline_app$yob
#calculate age (month) at censor (Apr0 2022)
baseline_app$diff_to_Apr <- (04 - baseline_app$mob)/12 
#add year and month together to get age at censor (Oct 2022) in years and months
baseline_app$Age_Apr_2022 <- baseline_app$diff_to_2022 + baseline_app$diff_to_Apr

#select columns of interest only - id, yob,mob,age at baseline, age oct 2020
baseline_app <- baseline_app[c(1,6,7,10,13)]

#############################################################################################################################################################################################################################

############################################################ CALCULATING AGE AT DEATH ####################################################################################################################################

#read in age at death file
age_dead <- read.csv("2023-02-01_dead.csv")

#check how many people died
dim(age_dead)
#[1] 1656    3

#filter to people who died prior to or in Apr 2022
age_dead_include <- age_dead[which(age_dead$ym <= 202204),]

#check how many people died in that timeframe
dim(age_dead_include)
#[1] 1656    3

#split date of death in year and month 
age_dead_include$y_dead <- as.numeric(substr(age_dead_include$ym, 1, 4))
age_dead_include$m_dead <- as.numeric(substr(age_dead_include$ym, 5, 6))

#merge date of death information with dob 
age_dead_include <- age_dead_include %>% left_join(dob, by = "id")

#calculate age (year) at death
age_dead_include$year_of_death <- age_dead_include$y_dead - age_dead_include$yob
#calculate age (month) at death
age_dead_include$month_of_death <- (age_dead_include$m_dead - age_dead_include$mob) / 12
#add year and month of death together to get age at death in years and months
age_dead_include$Age_at_death <- age_dead_include$year_of_death + age_dead_include$month_of_death 

#select columns of interest - id, age at death
age_dead_include <- age_dead_include[c(1,11)]

######################################################Creating death/alive files ######################################
#filter the baseline app file to indviduals who are alive
age_alive <- baseline_app[!baseline_app$id %in% age_dead_include$id,]

#check how many alive 
dim(age_alive)
#[1] 22435     3

#check the no. of people dead and alive add up to those in baseline app
22435 + 1656 
#[1] 24091

#it does match to baseline app 
dim(baseline_app)
#[1] 24091     3

#merge age at death information with baseline app information
age_dead_include <- age_dead_include %>% left_join(baseline_app, by = "id")

#create an age at censor column in death file set to NA for rbind purposes 
age_dead_include$Age_Apr_2022 <- NA

#create an age at death column in the alive file set to NA for rbind purposes 
age_alive$Age_at_death <- NA

#calculate time to event death in dead file by subtracting age at baseline from age at death 
age_dead_include$tte_death <- age_dead_include$Age_at_death - age_dead_include$Age_at_baseline

#calculate time to event death in alive file by subtracting age at baseline from age at censor 
age_alive$tte_death <- age_alive$Age_Apr_2022 - age_alive$Age_at_baseline

#create death event indicators in alive and dead files, dead = 1, alive = 0
age_dead_include$death_event <- 1
age_alive$death_event <- 0 

#reorganise columns in dead file to match order in alive file
age_dead_include <- age_dead_include[c(1, 3, 4, 5, 6, 2, 7, 8)]

#check names match 
names(age_alive)
#[1] "id"              "yob"             "mob"             "Age_at_baseline"
#[5] "Age_Apr_2022"    "Age_at_death"    "tte_death"       "death_event"
names(age_dead_include)
#[1] "id"              "yob"             "mob"             "Age_at_baseline"
#[5] "Age_Apr_2022"    "Age_at_death"    "tte_death"       "death_event"


#row bind the two datasets together 
all = rbind(age_alive, age_dead_include)

#check no. of rows in correct
dim(all)
#[1] 24091     8

#check how many dead and alive matches what we expect
table(all$death_event)

#0     1
#22435  1656

  

################################################ DEMENTIA ##############################################
#read in dementia diagnosis information 
DEM <- read.csv("all_dementia_prepped_15082023.csv")

#Check how many people have dementia 
dim(DEM)
#334 2

#split date of diagnosis information into month and year 
DEM$year_of_dementia = substring(DEM$first, 1, 4)
DEM$month_of_dementia = substring(DEM$first, 5,6)

#create dementia event indicator, dementia = 1 
DEM$dementia_event <- 1

#select columns of interest - id, year of dementia, month of dementia, dementia event 
age_onset = DEM[,c(1,3,4,5)]


#subset all file to those diagnosed with dementia to created an affected file
affected = all[which(all$id %in% DEM$id),] 

#check how many people in all file have dementia 
dim(affected)
#[1] 334   8  - everyone is present

#merge afftected file with age onset information 
affected = merge(age_onset, affected, by = "id")

#calculate age(month) at dementia diagnosis 
affected$month_diff_dem = (as.numeric(affected$month_of_dementia) - as.numeric(affected$mob))/12
#calculate age (year) at dementia diagnosis
affected$year_diff_dem = as.numeric(affected$year_of_dementia) - as.numeric(affected$yob)
# add together age year and month of diagnosis to get age at dementia diagnosis in years and months
affected$age_at_dementia = affected$year_diff_dem + affected$month_diff_dem

#calculate time to event for dementia by subtracting age at baseline from age at dementia diagnosis 
affected$tte_dementia <- affected$age_at_dementia - affected$Age_at_baseline
#filter age at dementia to 65 or over to avoid early onset dementia 
affected = affected[affected$age_at_dementia >=65,]

#how many cases were early onset 
dim(affected)
#[1] 312  15 

334 - 312 # 22 cases were early onset i.e. before age 65 


#check this is correct - confirmed these individuals were less than 65 years of age at dementia diagnosis
check <- age_onset[!age_onset$id %in% affected$id, ]
check <- check %>% left_join(dob, by ="id")
check$year <- as.numeric(check$year_of_dementia) - as.numeric(check$yob)
check$month <- as.numeric(check$month_of_dementia) - as.numeric(check$mob)/12


#remane id to sample name
colnames(affected)[1] = "Sample_Name"

#select columns of interest only - sample name, dementia event, tte dementia, death event, tte death, age at baseline 
affected_final <- affected[c(1,4,15,11,10,7)]

#check for missing info 
table(is.na(affected_final))
#FALSE
#1872

#subset all file to indviduals who were not diagnosed with dementia on or before Apr 2022
healthy = all[-which(all$id %in% DEM$id),]

#create dementia event, no dementia = 0
healthy$dementia_event <- 0

#calculate time to event dementia for dementia free individuals 
healthy$tte_dementia <- ifelse(healthy$death_event == 1, (healthy$Age_at_death - healthy$Age_at_baseline), (healthy$Age_Apr_2022 - healthy$Age_at_baseline))

#create an age filter to create an age appropriate control group of non dementia partcipants 
healthy$age_filter <- ifelse(healthy$death_event == 1, healthy$Age_at_death, healthy$Age_Apr_2022)

#use age filter column to remove indviduals who were younger than 65 in Oct 2020 or died before age 65 
healthy = healthy[healthy$age_filter >=65,]

#rename id to sample name
colnames(healthy)[1] <- "Sample_Name"

#select columns of interest - sample name, dementia event, tte dementia, death event, tte death, age at baseline 
healthy_final <- healthy[c(1, 9, 10, 8, 7, 4)]

#check for missing info - 11 NA's for the indviduals with missing age and sex data as mentioned above
summary(healthy_final)

#remove individuals with missing age and sex information
healthy_final <- healthy_final[complete.cases(healthy_final), ]

#check for NA's 
table(is.na(healthy_final))
#FALSE
#59328


#check column names match
names(affected_final)
#[1] "Sample_Name"     "dementia_event"  "tte_dementia"    "death_event"
#[5] "tte_death"       "Age_at_baseline"
names(healthy_final)
#[1] "Sample_Name"     "dementia_event"  "tte_dementia"    "death_event"
#[5] "tte_death"       "Age_at_baseline"


#row bind affected and healthy datasets together 
combined <- rbind(affected_final, healthy_final)

#check ids are all unique 
length(unique(combined$Sample_Name))
#[1] 10200

dim(combined)
#[1] 10200     6

#read in GS target file with methylation info
target <- readRDS("GS20k_Targets.rds")
#subset to indivduals with methylation data
combined <- combined[which(combined$Sample_Name %in% target$Sample_Name),] 

#check how many with methylation
dim(combined)
#[1] 7791    6

#read in EpiScores
prot <- read.csv("KORA_Episcores_projected_adjusted_in_GS20K_19012023.csv")
#Select columns of interest
prot <- prot[c(1:86)]

#combined time to event data with episcores
combined_epi <- combined %>% left_join(prot, by = "Sample_Name")
#select some of the target info
target <- target[c(1,4)]
#combined columns of interest from target file with time to event data
combined_epi <- combined_epi %>%  left_join(target, by = "Sample_Name")

#recode sex to binary rather than characters 
combined_epi <- mutate(combined_epi, sex=recode(sex, 'M' = '0', 'F' = '1'))

#save of time to event out
write.csv(combined_epi, file = "time_to_event_data_dementia_death_GenScot_84_EpiScores_15082023.csv")

