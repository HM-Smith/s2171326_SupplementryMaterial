#read in data 
library(dplyr)
data <- read.csv("LBC1936_WAVE1_cohort_data_KORA_EpiScore_Adjusted_17072023.csv")

#set sex as factor 
data$sex <- as.factor(data$sex)

#select columns of interest
data <- data %>% select(3:65, 73, 148, 174:188,213:229, 328:332, 335:419, 153 )

#read in and combine with epismoker 
epismoke <- readRDS("lbc_epismoker.rds") 
names(epismoke)[1] <- "Basename"
data <- data %>% left_join(epismoke, by = "Basename")

#99999 = missing code
summary(data$depind_w1)
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
#4    3134    5354    5265    6264   99999

#recode missing code to NA 
data$depind_w1 <- ifelse(data$depind_w1 == 99999, NA, data$depind_w1)

#check recode worked
summary(data$depind_w1)
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's
# 4    3092    5342    4565    6258    6505       8

#select cognitive tests 
covariates <- colnames(data)[c(11:14, 68, 15:18, 69, 19:22, 70, 23:26, 71, 27: 30, 72, 31:34, 73, 35:38, 74, 39:42, 75, 43:46, 76, 47:50, 77, 51:54, 78, 55:58, 79, 59:62, 80)]

#create empty dataframe for stats
desc_stats <- data.frame(Test = covariates, N = NA, Mean = NA, SD = NA, Missing = NA, Range =  NA)

#set cog test names as rownames
rownames(desc_stats) <- covariates

#loop for looping over cognitive tests and calculating descriptive stats
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
  desc_stats[i, 6] <- paste("'",round(Range[1], 2), "-", round(Range[2],2))
  
}

#create cohort, test and wave columns 
desc_stats$Cohort <- rep("LBC1936", 65)
desc_stats$Test <- c("Digit symbol backwards", "Digit symbol backwards", "Digit symbol backwards", "Digit symbol backwards","Digit symbol backwards",
                          "Matrix Reasoning","Matrix Reasoning","Matrix Reasoning","Matrix Reasoning","Matrix Reasoning",
                          "Block design", "Block design","Block design","Block design","Block design",
                          "Digit symbol", "Digit symbol","Digit symbol","Digit symbol","Digit symbol",
                          "Symbol search",  "Symbol search",  "Symbol search",  "Symbol search",  "Symbol search", 
                          "National adult reading test", "National adult reading test", "National adult reading test", "National adult reading test", "National adult reading test", 
                          "Logical memory","Logical memory","Logical memory","Logical memory","Logical memory",
                          "Verbal paired associates", "Verbal paired associates", "Verbal paired associates", "Verbal paired associates", "Verbal paired associates", 
                          "Inspection time", "Inspection time","Inspection time","Inspection time","Inspection time",
                          "Four choice reaction time", "Four choice reaction time","Four choice reaction time","Four choice reaction time","Four choice reaction time",
                          "Spatial span", "Spatial span","Spatial span","Spatial span","Spatial span",
                          "Verbal fluency",  "Verbal fluency", "Verbal fluency", "Verbal fluency", "Verbal fluency",
                          "Wechsler test of reading","Wechsler test of reading","Wechsler test of reading","Wechsler test of reading","Wechsler test of reading" 
                          )

desc_stats$Wave <- rep(1:5, 13)

#save out files
write.csv(desc_stats, file = "LBC1936_Cog_Tests_Desc_Stats_07082023.csV", row.names = F)

#########################################################################################################################################################################################################
## select covariates
 covariates <- colnames(data)[c(98:102,63:65, 188:189)]

#create empty dataframe for stats
desc_stats <- data.frame(Covariate = covariates, N = NA, Mean = NA, SD = NA, Missing = NA, Range = NA)

#set covariate names as rownames 
rownames(desc_stats) <- covariates

#loop for looping over covariate tests and calculating descriptive stats
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
  desc_stats[i, 6] <- paste("'",round(Range[1], 2), "-", round(Range[2],2))
  
}

#create cohort, covariate and wave columns 
desc_stats$Cohort <- rep("LBC1936", 10)
desc_stats$Covariate <- c("Age", "Age", "Age", "Age", "Age",
                          "Body Mass Index", 
                          "Alcohol units (per week)", 
                          "Scottish deprivation index rank",
                          "Years in full time education", 
                          "Epigenetic smoking score")
desc_stats$Wave <- c(1:5, 1, 1, 1, 1, 1)


#save out file
write.csv(desc_stats, file = "LBC1936_Covariates_Tests_Desc_Stats_27102023.csV", row.names = F)

#number of males and females 
table(data$sex) # 0 = male, 1 = female 
#0   1
#548 543

#percentage of males and females 
548/(548+543)*100
#50.22915

543/(548+543)*100
#[1] 49.77085

#########################################################################################################################################################################################################
## MRI 
#create list of MRI variable names
covariates <- colnames(data)[c(81:97)]

#create empty dataframe for stats
desc_stats <- data.frame(Measure = covariates, N = NA, Mean = NA, SD = NA, Missing = NA, Range = NA)

#set MRI variable names as rownames
rownames(desc_stats) <- covariates


#loop for looping over MRI variables and calculating descriptive stats
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
  desc_stats[i, 6] <- paste("'",round(Range[1], 2), "-", round(Range[2],2))
}


#create cohort, measure and wave columns 
desc_stats$Cohort <- rep("LBC1936", 17)
desc_stats$Measure <- c("Intrcranial volume",
                      "Total brain volume","Total brain volume","Total brain volume","Total brain volume",
                      "White matter hyperintensity volume", "White matter hyperintensity volume","White matter hyperintensity volume","White matter hyperintensity volume",
                      "Grey matter volume", "Grey matter volume","Grey matter volume","Grey matter volume",
                      "Normal appearing white matter volume","Normal appearing white matter volume","Normal appearing white matter volume","Normal appearing white matter volume")

desc_stats$Wave <- c(2, rep(2:5, 4))

#save out file
write.csv(desc_stats, file = "LBC1936_MRI_Measures_Desc_Stats_07082023.csV", row.names = F)
