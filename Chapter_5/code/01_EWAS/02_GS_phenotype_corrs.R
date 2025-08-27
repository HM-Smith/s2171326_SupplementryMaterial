#############################################################

## Investigating traits and demographics 

#############################################################

## load libraries
library(tidyverse)

library(dplyr)
library(coxme)
library(kinship2)
library(data.table)


## MAKING DATASET

## read in body data 
body <- read.csv("body.csv")

names(body)[1] <- "ID"

## view subsection of data 
head(body)

## read in covariates
covars <- read.csv("covariates.csv")

names(covars)[1] <- "ID"

covars <- covars %>% select(ID, age, sex, Sample_Sentrix_ID)

## check covars 
head(covars)

## read in biochemistry data
biochem <- read.csv("biochemistry.csv")

## check subsection of biochemistry data 
head(biochem)



## merge datasets by ID 
combined_data <- covars %>%  left_join(body, by = "ID")
combined_data <- combined_data %>% left_join(biochem, by = "ID")



ped1 <- read.csv("/Cluster_Filespace/Marioni_Group/Hannah/Biomarkers_brain_health_LBC_GS/GenScot_input_data/pedigree_data/2023-03-20_pedigree.csv")
ped2 = data.frame(famid=c(4091,4384), volid=c(103027, 144865), father=c(0,0), mother=c(0,0), sex= c("F", "F"))
ped = rbind(ped1, ped2)
kin <- with(ped, pedigree(volid, father, mother, sex, famid=famid))
kin_model <- kinship(kin) 



## Subset cov file to those in pedigree file to avoid error in lmekin function (can't have anyone in phenotype file thats not in kinship matrix)
cov = combined_data[which(combined_data$ID %in% ped$volid),]


dim(combined_data)
#[1] 18413    20


dim(cov)
#[1] 18411    20


#####  filtering the cov file for Ids that are present in the residualised phenotype files ####
bmi_res <- read.csv("bmi_res.csv")

bmi_include <- cov[cov$Sample_Sentrix_ID %in% bmi_res$ID,]
bmi_include <- bmi_include %>% select(Sample_Sentrix_ID, ID, bmi)


WHR_res <- read.csv("whr_res.csv")

WHR_include <- cov[cov$Sample_Sentrix_ID %in% WHR_res$ID,]
WHR_include <- WHR_include %>% select(ID, whr)


body_res <- read.csv("body_fat_res.csv")

body_include <- cov[cov$Sample_Sentrix_ID %in% body_res$ID,]
body_include <- body_include %>% select(ID, body_fat)


glucose_res <- read.csv("Glucose_res.csv")

glucose_include <- cov[cov$Sample_Sentrix_ID %in% glucose_res$ID, ]
glucose_include <- glucose_include %>% select(ID, Glucose)

HDL_res <- read.csv("HDL_cholesterol_res.csv")

HDL_include <- cov[cov$Sample_Sentrix_ID %in% HDL_res$ID, ]
HDL_include <- HDL_include %>% select(ID, HDL_cholesterol)


total_res <- read.csv("Total_cholesterol_res.csv")

total_include <- cov[cov$Sample_Sentrix_ID %in% total_res$ID,]
total_include <- total_include %>% select(ID, Total_cholesterol)

##########################################################################################

#select age and id
age <- cov %>% select(ID, age)

#merge datasets together
data <- full_join(bmi_include, WHR_include, by = "ID")
data <- full_join(data, body_include, by = "ID")
data <- full_join(data, glucose_include, by = "ID")
data <- full_join(data, HDL_include, by = "ID")
data <- full_join(data, total_include, by = "ID")
data <- full_join(data, age, by = "ID")

#read in meth pcs
quant <- read.table("methPC20_quant_18411.cov", header = T, sep=' ')
names(quant)[1] <- "Sample_Sentrix_ID"

#meth meth pcs into dataset
data <- full_join(data, quant, by = "Sample_Sentrix_ID")
data <- data %>% select(-IID, -ID, -Sample_Sentrix_ID) 

#rename covariates and phenotypes for plot
names(data) <- c("BMI (kg/m2)", "WHR", "Body fat (%)", "Glucose (mmol/L)", "HDL cholesterol (mmol/L)", "Total cholesterol (mmol/L)", "age (years)", "EpiSmokEr", "CD8+ T cells",  "CD4+ T cells",
                 "Natural Killer cells", "Monocytes", "B cells", "PC1", "PC2", "PC3", "PC4", "PC5",  "PC6", "PC7", "PC8","PC9", "PC10", "PC11", "PC12", "PC13", "PC14",  "PC15", "PC16", "PC17", "PC18",
                 "PC19","PC20")

#select appropriate columns
pheno_cov <- data[, 1:13]

#make correlation matrix of covars and phenotypes
res <- cor(pheno_cov, use="pairwise.complete.obs")

#save out correlation matrix
write.csv(res, "GS_phenotype_covar_corr.csv")

#load corrplot library
library(corrplot)

#plot correlation matrix
setwd("")
tiff("GS_phenotype_covar_corr.tiff", width = 8*300, height = 9*300, res = 300)

#par(mar = c(2, 2, 2, 2))

corrplot(res, type = "lower", order = "hclust", tl.col = "black", addCoef.col = 'black', method = "color", number.digits = 2)
dev.off()

#make correlation matrix of data and set as dataframe
res <- cor(data, use="pairwise.complete.obs")
df <- as.data.frame(res)

#select Pc columns and set a maxtrix 
df2 <- df[1:13,14:ncol(df)]
df2 <- as.matrix(df2)

#plot the correlation matrix of PCs
setwd("")
tiff("GS_phenotype_covar_DNAmPCscorr.tiff", width = 14*300, height = 10*300, res = 300)

#par(mar = c(2, 2, 2, 2))

corrplot(df2, tl.col = "black", addCoef.col = 'black', method = "color", number.digits = 2)
dev.off()


