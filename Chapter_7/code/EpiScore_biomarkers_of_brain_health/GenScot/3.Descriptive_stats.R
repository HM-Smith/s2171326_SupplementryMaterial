library(dplyr)

#read in GS data 
data <- read.csv("GS20K_Cognitive_KORA_adjust_EpiScores_combined_20012023.csv")

#set sex as factor 
data$sex <- as.factor(data$sex)

#select cogntive tests
covariates <- colnames(data)[c(6,10:11, 23)]

#create empty dataframe 
desc_stats <- data.frame(Covariate = covariates, N = NA, Mean = NA, SD = NA, Missing = NA, Range = NA)

#set cognitive tests as rownames
rownames(desc_stats) <- covariates

#loop over cognitive tests and calculate descriptive statistics 
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
  desc_stats[i, 6] <- paste0("'",Range[1], "-", Range[2])
  }

#create cohort and covariate columns 
desc_stats$Cohort <- rep("Generation Scotland", 4)
desc_stats$Covariate <- c("Digit symbol", "Verbal total", "Vocabulary", "Logical memory")

#save out file
write.csv(desc_stats, file = "GenScot_Cog_Tests_Desc_Stats_07082023.csv", row.names = F)

#########################################################################################################################################################################################################
library(tidyverse)

## Covariates 

#read in GS data
data <- read.csv("GS20K_Cognitive_KORA_adjust_EpiScores_combined_20012023.csv")

#set sex as factor 
data$sex <- as.factor(data$sex)

#read in covariate files
alc <- read.csv("2023-02-01_alc.csv")
smid <- read.csv("2023-02-13_simd.csv")
BMI <- read.csv("body.csv")

#rename id column 
names(alc)[1] <- "Sample_Name"
names(smid)[1] <- "Sample_Name"
names(BMI)[1] <- "Sample_Name"

#read in EpiSmokEr scores and combine
w1 <- readRDS("wave1_epismoker.rds")
w3 <- readRDS("wave3_epismoker.rds")
w4 <- readRDS("wave4_epismoker.rds")
bind <- rbind(w1, w3)
bind <- rbind(bind, w4) # 18779 individuals with DNAm epismoker calculated 
bind$Sample_Sentrix_ID <- row.names(bind)
bind <- bind[-1]

#join cognitive data and covariates
data <- data %>% left_join(alc, by ="Sample_Name")
data <- data %>% left_join(smid, by ="Sample_Name")
data <- data %>% left_join(BMI, by ="Sample_Name")
data <- data %>% left_join(bind, by ="Sample_Sentrix_ID")

#select columns of interest
data <- data[, c(109,113,115,117,120,125,110)]

#select covariates
covariates <- colnames(data)[c(1:6)]

#create empty dataframe
desc_stats <- data.frame(Covariate = covariates, N = NA, Mean = NA, SD = NA, Missing = NA, Range = NA)

#set covariates as rownames
rownames(desc_stats) <- covariates

#for loop that calculates descriptive statistics 
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

#create cohort and covariate columns 
desc_stats$Cohort <- rep("Generation Scotland", 6)
desc_stats$Covariate <- c("Age", "Alcohol units (per week)", "Years in full time education", "Scottish deprivation index rank", "Body Mass Index", "Epigenetic smoking score")

#save out file
write.csv(desc_stats, file = "GenScot_Covariates_Tests_Desc_Stats_07082023.csv", row.names = F)


#check number of males and females
table(data$sex) # M= 1, F= 2
#1     2
#7580 10833


#calculate percentage of males and females
7580 / (7580+10833)*100

#[1] 41.16657

10833 / (7580+10833)*100
#[1] 58.83343


