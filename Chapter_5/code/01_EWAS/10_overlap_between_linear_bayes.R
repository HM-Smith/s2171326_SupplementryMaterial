library(dplyr)

#read in OSCA significant results files
BMI <- read.csv("bmi_EWAS_20methPCs_18411_significant.csv")
Body <- read.csv("body_fat_EWAS_20methPCs_18411_significant.csv")
WHR <- read.csv("whr_EWAS_20methPCs_18411_significant.csv")
Glucose <- read.csv("Glucose_EWAS_20methPCs_18411_significant.csv")
HDL <- read.csv("HDL_cholesterol_EWAS_20methPCs_18411_significant.csv")
Total <- read.csv("Total_cholesterol_20methPCs_18411_significant.csv")


#read in significant bayesian results
bayes <- read.csv("Bayes_Results/significant_hits_full_model.csv")

#check number of hits per trait in bayesian file
table(bayes$trait)
#bmi          body_fat           Glucose   HDL_cholesterol
#27                18                 3                20
#Total_cholesterol               whr
#19                12

#filter to BMI results
bmi_common <- bayes[bayes$trait == "bmi",]

#check num of hits
dim(bmi_common)
#[1] 27  9

#check overlap with OSCA bmi results
bmi_common <- bayes[bayes$trait == "bmi" & bayes$Name %in% BMI$Probe,]

#check num of common hits
dim(bmi_common)
#[1] 25  9

#filter to body fat results 
body_fat_common <- bayes[bayes$trait == "body_fat",]

#check num of sig hits
dim(body_fat_common)
#[1] 18  9

#check overlap with OSCA 
body_fat_common <- bayes[bayes$trait == "body_fat" & bayes$Name %in% Body$Probe,]

#check num of common hits
dim(body_fat_common)
#[1] 17  9


#filter to glucose sig results
Glucose_common <- bayes[bayes$trait == "Glucose",]

#check num of sig hits
dim(Glucose_common)
#[1] 3  9

#check overlap with OSCA results
Glucose_common <- bayes[bayes$trait == "Glucose" & bayes$Name %in% Glucose$Probe,]

#check num of overlapping cpgs
dim(Glucose_common)
#[1] 1 9


#filter to whr sig results
WHR_common <- bayes[bayes$trait == "whr",]

#check num of sig hits
dim(WHR_common)
#[1] 12  9

#check overlap with OSCA results
WHR_common <- bayes[bayes$trait == "whr" & bayes$Name %in% WHR$Probe,]

#check num of overlapping cpgs
dim(WHR_common)
#[1] 11  9

 
#filter to HDL cholesterol sig hits 
HDL_common <- bayes[bayes$trait == "HDL_cholesterol",]

#check num of sig hits
dim(HDL_common)
#[1] 20  9

#check overlap with OSCA results
HDL_common <- bayes[bayes$trait == "HDL_cholesterol" & bayes$Name %in% HDL$Probe,]

#check number of overlapping cpgs
dim(HDL_common)
#[1] 16  9


#filter to total cholesterol results
total_common <- bayes[bayes$trait == "Total_cholesterol",]

#check num of sig hits
dim(total_common)
#[1] 19  9

#check num of overlap with OSCA results
total_common <- bayes[bayes$trait == "Total_cholesterol" & bayes$Name %in% Total$Probe,]

#check num of overlapping cpgs
dim(total_common)
#[1] 18  9


