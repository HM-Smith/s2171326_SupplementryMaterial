library(dplyr)

#read in LBC1921 data
data <- read.csv("LBC1921_WAVE1_cohort_data_KORA_EpiScore_Adjusted_31012023.csv")

#set sex as factor
data$sex <- as.factor(data$sex)

#select columns of interest
data <- data %>% select(3, 9:27, 32, 41:42, 49, 102:189, 193)

#read in epismoker and rename id column 
epismoke <- readRDS("lbc_epismoker.rds") 
names(epismoke)[1] <- "Basename"

#merge data and episcore 
data <- data %>% left_join(epismoke, by = "Basename")


#select cognitive test names
covariates <- colnames(data)[c(2:20)]

#create empty dataframe 
desc_stats <- data.frame(Test = covariates, N = NA, Mean = NA, SD = NA, Missing = NA, Range = NA)

#set cogntive test names as rownames
rownames(desc_stats) <- covariates

#loop over cognitive tests calculating descriptive statistics 
for(i in covariates){
  
  m <- mean(data[,i], na.rm = T)
  SD <- sd(data[,i], na.rm =T)
  missing <- sum(is.na(data[,i]))
  n <- sum(!is.na(data[,i]))
  Range <- range(data[,i], na.rm = T)
  
  desc_stats[i, 1] <- as.character(i)
  desc_stats[i, 2] <- n
  desc_stats[i, 3] <- m
  desc_stats[i, 4] <- SD
  desc_stats[i, 5] <- missing
  desc_stats[i, 6] <- paste0("'",round(Range[1], 2), "-", round(Range[2],2))
  
}

#rewrite names of cognitive tests
desc_stats$Test <- c("Verbal fluency",  "Verbal fluency", "Verbal fluency", "Verbal fluency", "Verbal fluency",
                     "Ravens matrices", "Ravens matrices", "Ravens matrices", "Ravens matrices", "Ravens matrices",
                     "National adult reading test", "National adult reading test", "National adult reading test", "National adult reading test", 
                     "Logical memory","Logical memory","Logical memory","Logical memory","Logical memory")

#set waves and cohort indicator 
desc_stats$Wave <- c(1:5, 1:5, 1:4, 1:5)
desc_stats$Cohort <- rep("LBC1921", 19)

#select columns of interest
desc_stats <- desc_stats %>% select(1,6,2:5,7)

#save out cognitive test descriptive stats 
write.csv(desc_stats, file = "LBC1921_Cog_Tests_Desc_Stats_07082023.csv", row.names = F)
#########################################################################################################################################################################################################
##Covariates

#select covariate names
covariates <- colnames(data)[c(113, 21:24,114)]

#create empty dataframe 
desc_stats <- data.frame(Covariate = covariates, N = NA, Mean = NA, SD = NA, Missing = NA, Range = NA)

#set covariate names as rownames
rownames(desc_stats) <- covariates


#loop over covariates and calculate descriptive stats 
for(i in covariates){n
  
  m <- mean(data[,i], na.rm = T)
  SD <- sd(data[,i], na.rm =T)
  missing <- sum(is.na(data[,i]))
  n <- sum(!is.na(data[,i]))
  Range <- range(data[,i], na.rm = T)
  
  desc_stats[i, 1] <- as.character(i)
  desc_stats[i, 2] <- n
  desc_stats[i, 3] <- m
  desc_stats[i, 4] <- SD
  desc_stats[i, 5] <- missing
  desc_stats[i, 6] <- paste0("'", round(Range[1],2), "-", round(Range[2], 2))
  
}

#rewrite names of covariates 
desc_stats$Covariate <- c("Age", 
                          "Alcohol units (per week)", 
                          "Scottish deprivation index rank",
                          "Years in full time education", 
                          "Body Mass Index", 
                          "Epigenetic smoking score")

#create wave and cohort indicators 
desc_stats$Wave <- rep(1, 6)
desc_stats$Cohort <- rep("LBC1921", 6)
#desc_stats <- desc_stats %>% select(1,6, 2:5, 7)

#save out descriptive stats for covariates 
write.csv(desc_stats, file = "LBC1921_Covariates_Tests_Desc_Stats_07082023.csV", row.names = F)

#check number of males and females
table(data$sex) # 0 = male, 1 = female 
#0   1
#238 331

#calculate percentages 
238/(238+331)*100
#[1] 41.82777

331/(238+331)*100
#[1] 58.17223


#########################################################################################################################################################################################################

#check dementia cases 
table(data$dement_consensus)

#       M   N   Y
# 439   7  13 110

table(data$dement_consensus, data$dement_subtype)

#         MCI POSSIBLE MCI POSSIBLE PDD POSSIBLE VD PROBABLE AD PROBABLE MIX
#  439   0            0            0           0           0            0
#M   0   0            0            0           1           0            0
#N  10   2            1            0           0           0            0
#Y   0   0            0            1           6          38            9

#            PROBABLE PSP PROBABLE VD UNKNOWN
#             0           0       0
#M            0           0       6
#N            0           0       0
#Y            1          25      30


table(data$dement_consensus, data$sex)

#0   1
#187 252
#M   3   4
#N   7   6
#Y  41  69
