#read in cognitive tests 
cog <- read.csv("cognitive_tests_GS_08102024.csv")

#check data dimensions 
dim(cog)
#[1] 21510     5

#check column names
names(cog)
#[1] "Sample_Name"  "digit_symbol" "verbal_total" "vocabulary"   "LM"

#rename sample name to id 
names(cog)[1] <- "id"

#load data.table package
library(data.table)

#read in AD PRS scores 
AD_PRS <- fread("alzheimers_aggregated_scores.txt.gz")

#select ID and scores 
AD_PRS <- AD_PRS[, c(2, 4)]

#rename id column for merging 
names(AD_PRS)[1] <- "id"

#load dplyr 
library(dplyr)

#join AD PRS with cog tests 
df <- full_join(cog,AD_PRS, by = "id")

#read in cognitive PRS
cog_PRS <- fread("aggregated_scores.txt.gz")

#select id and scores 
cog_PRS <- cog_PRS[,c(2,4)]

#rename id column for merging 
names(cog_PRS)[1] <- "id"

#merge with cog tests and AD PRS
df<- full_join(df,cog_PRS, by = "id")

#remove Na's - as complete data is required for MAJA 
df <- df[complete.cases(df), ]

#check dimensions 
dim(df)
#[1] 19313     7

#read in protein data 
prot <- readRDS("GS_ProteinGroups_RankTransformed_23Aug2023.rds")

#check prot dimensions 
dim(prot)
#[1] 15818   440

#join cog and proteins 
df <- inner_join(df, prot, by = "id")

#check dimensions
dim(df)
#[1] 14537   446


#read in covars 
covars <- read.csv("covars_for_cog_post_outlier_removal_KNNimputed.csv")

#check column names of covarsfile 
names(covars)
# [1] "id"           "rank"         "bmi"          "units"        "age"
# [6] "sex"          "depression_Y" "years"        "diabetes_Y"   "high_BP_Y"
# [11] "smokingScore" "e4_count"

dim(covars)
#[1] 14258    13

#join in the covars
df <- inner_join(df, covars, by = "id")


#check data dimensions 
dim(df)
#[1] 14258   468


#remove Na's - as complete data is required for MAJA 
df <- df[complete.cases(df), ]

#check data dimensions
dim(df)
#[1] 14258   468

write.csv(df, file= "complete_data_KNN_imputed.csv" )

#check NA's 
table(is.na(df))
# FALSE
# 6530164

#read in coxme and kinship packages 
library(coxme)
library(kinship2)

#read in pedigree file 
ped <- read.csv("2023-03-20_pedigree.csv")

#make kinship matrix 
kin <- with(ped, pedigree(volid, father, mother, sex, famid=famid))
kin_model <- kinship(kin) 


#save list of protein names
protein_list <- colnames(df)[8:446]

#check length of protein list
length(protein_list)
#[1] 439


#load foreach and doParallel packages 
library(foreach)
library(doParallel)

# detect and register cores
cores <- detectCores()
cl <- makeCluster(20)
registerDoParallel(cl)

# parallelised for each loop that adjusts protein levels for covars, saves out dataframe of residuals 
x <- foreach(i = protein_list, .packages = c("coxme"), .combine='cbind') %dopar% { 
  mod <- lmekin(df[, i] ~ scale(age)  + factor(sex)  + scale(units) + scale(bmi) + scale(rank) + scale(smokingScore) + factor(depression_Y) + 
                  scale(years) + factor(high_BP_Y) + factor(diabetes_Y) + factor(e4_count) + (1|id), varlist = kin_model*2, data = df, na.action = na.exclude)
  protein_res <- resid(mod)
  protein_res <- as.data.frame(protein_res)
  names(protein_res)  <- i 
  result <- protein_res
}

#check dimensions of x (residuals from adjusted model)
dim(x)
#[1] 14258   439


#stop cluster cores
stopCluster(cl)


#scale residuals 
for(protein in protein_list){
  y <- scale(x[,protein])
  x[,protein] <- as.vector(y) 
}


#load readr package
library(readr)

#write out age, sex and relatedness corrected and scaled proteins 
write_rds(x, "GS_MS_prots_MAJA_cog_adPRS_cogPRS_KNNIMPUTE_corrected.rds")

#select outcome columns 
cog <- df[,2:7]

#print column names 
names(cog)
# [1] "digit_symbol"             "verbal_total"
# [3] "vocabulary"               "LM"
# [5] "alzheimers_scorefile_SUM" "cognition_scorefile_SUM"



#check cog dimensions
dim(cog)
#[1] 14258     6
  

#scale outcomes 
cog$digit_symbol <- scale(cog$digit_symbol)
cog$verbal_total <- scale(cog$verbal_total)
cog$vocabulary <- scale(cog$vocabulary)
cog$LM <- scale(cog$LM)
cog$alzheimers_scorefile_SUM <- scale(cog$alzheimers_scorefile_SUM)
cog$cognition_scorefile_SUM <- scale(cog$cognition_scorefile_SUM)


#save out outcomes to text file
write.table(cog, "GS_cognitive_adPRS_cogPRS_MAJA_KNNimpute_corrected.txt", row.names = FALSE, col.names = FALSE)


#############################################################################################################################################################################################################

# convert protein data to .zarr file with python
eval "$(micromamba shell hook --shell bash )"
micromamba activate MAJA
source maja-python/bin/activate maja-python

# #Install required packages
# pip install -U pip
# pip install numpy scipy matplotlib loguru mpi4py welford zarr dask pandas tqdm bed_reader scikit-allel pyreadr

python 

import numpy as np
import zarr
import pyreadr

df = pyreadr.read_r("GS_MS_prots_MAJA_cog_adPRS_cogPRS_KNNIMPUTE_corrected.rds")

df = df[None]

df.shape
(14258, 439)



arr = df.to_numpy()

arr.shape
(14258, 439)


z = zarr.array(arr, store = "GS_MS_prots_MAJA_cog_adPRS_cogPRS_KNNIMPUTE_corrected.zarr")