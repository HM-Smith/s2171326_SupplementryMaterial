#read in cognitive tests 
cog <- read.csv("")

#check data dimensions 
dim(cog)
#[1] 21510     5

#check column names
names(cog)
#[1] "Sample_Name"  "digit_symbol" "verbal_total" "vocabulary"   "LM"

#rename sample name to id 
names(cog)[1] <- "id"

#read in age and sex info 
age_sex <- read.csv("")

#check column names of age and sex file 
names(age_sex)
#[1] "id"  "age" "sex"

#load dplyr library
library(dplyr)

#merge together cog tests and age & sex info 
df <- left_join(cog, age_sex, by = "id")


#check for NA's
table(is.na(df$age))
#FALSE
#21510

table(is.na(df$sex))
#FALSE
#21510

table(is.na(df$digit_symbol))
#FALSE  TRUE
#21220   290

table(is.na(df$LM))
#FALSE  TRUE
#21207   303

table(is.na(df$verbal_total))
#FALSE  TRUE
#21200   310


table(is.na(df$vocabulary))
#FALSE  TRUE
#21038   472

#remove Na's - as complete data is required for MAJA 
df <- df[complete.cases(df), ]

#check NA's 
table(is.na(df))
#FALSE
#144998

dim(df)
#[1] 20714     7

#load data.table package
library(data.table)

#read in AD PRS scores 
AD_PRS <- fread("")

#select ID and scores 
AD_PRS <- AD_PRS[, c(2, 4)]

#rename id column for merging 
names(AD_PRS)[1] <- "id"

#join AD PRS with cog tests 
df <- df %>% left_join(AD_PRS, by = "id")

#read in cognitive PRS
cog_PRS <- fread("")

#select id and scores 
cog_PRS <- cog_PRS[,c(2,4)]

#rename id column for merging 
names(cog_PRS)[1] <- "id"

#merge with cog tests and AD PRS
df<- df %>% left_join(cog_PRS, by = "id")

#remove Na's - as complete data is required for MAJA 
df <- df[complete.cases(df), ]

#check for NA's 
table(is.na(df))
#FALSE
#173817

#check df dimensions 
dim(df)
#[1] 19313     9

#read in protein data 
prot <- readRDS("")

#check prot dimensions 
dim(prot)
#[1] 15818   440

#merge prot with cog tests and PRSs
df <- prot %>% left_join(df, by = "id")

#check for NA's 
table(is.na(df))
#FALSE    TRUE
#7076216   10248

#remove Na's - as complete data is required for MAJA 
df <- df[complete.cases(df), ]

#check for NA's 
table(is.na(df))
#FALSE
#6512576

#check df dimensions
dim(df)
#[1] 14537   448

#save out complete dataset 
write.csv(df, "")

#read in coxme and kinship packages 
library(coxme)
library(kinship2)

#read in pedigree file 
ped <- read.csv("")

#make kinship matrix 
kin <- with(ped, pedigree(volid, father, mother, sex, famid=famid))
kin_model <- kinship(kin) 


#save list of protein names
protein_list <- colnames(df)[2:440]

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

# parallelised for each loop that adjusts protein levels for age, sex and relatedness, saves out dataframe of residuals 
x <- foreach(i = protein_list, .packages = c("coxme"), .combine='cbind') %dopar% { 
  mod <- lmekin(df[, i] ~ scale(age)  + factor(sex)  + (1|id), varlist = kin_model*2, data = df, na.action = na.exclude)
  protein_res <- resid(mod)
  protein_res <- as.data.frame(protein_res)
  names(protein_res)  <- i 
  result <- protein_res
}

#check dimensions of x (residuals from adjusted model)
dim(x)
#[1] 14537   439

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
write_rds(x, "")

#select outcome columns 
cog <- df[,c(441:444, 447:448)]

#print column names 
names(cog)
#[1] "digit_symbol"             "verbal_total"
#[3] "vocabulary"               "LM"
#[5] "alzheimers_scorefile_SUM" "cognition_scorefile_SUM"

#check cog dimensions
dim(cog)
#[1] 14537     6


#scale outcomes 
cog$digit_symbol <- scale(cog$digit_symbol)
cog$verbal_total <- scale(cog$verbal_total)
cog$vocabulary <- scale(cog$vocabulary)
cog$LM <- scale(cog$LM)
cog$alzheimers_scorefile_SUM <- scale(cog$alzheimers_scorefile_SUM)
cog$cognition_scorefile_SUM <- scale(cog$cognition_scorefile_SUM)


#save out outcomes to text file
write.table(cog, "", row.names = FALSE, col.names = FALSE)


#############################################################################################################################################################################################################

##covert proteins to .zarr files in python 
eval "$(micromamba shell hook --shell bash )"
micromamba activate BayesR
source BayesR_python/bin/activate


python 

import numpy as np
import zarr
import pyreadr

df = pyreadr.read_r("")

df = df[None]

df.shape
(14537, 439)


arr = df.to_numpy()

arr.shape
(14537, 439)


z = zarr.array(arr, store = "")



