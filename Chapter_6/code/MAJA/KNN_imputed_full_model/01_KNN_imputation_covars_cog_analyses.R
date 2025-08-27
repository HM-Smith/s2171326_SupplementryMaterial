#load dplyr 
library(dplyr)
#read in cognitive test data 
cog <- read.csv("cognitive_tests_GS_08102024.csv")

#subset to complete cases 
cog <- cog[complete.cases(cog),]

#load data table library 
library(data.table)

#read in cognitvie PRS and select columns 
cog_PRS <- fread("aggregated_scores.txt.gz")
cog_PRS <- cog_PRS[,c(2,4)]
names(cog_PRS)[1] <- "Sample_Name"

#read in AD_PRS and select columns 
AD_PRS <- fread("alzheimers_aggregated_scores.txt.gz")
AD_PRS <- AD_PRS[, c(2, 4)]
names(AD_PRS)[1] <- "Sample_Name"

#join cognitive tests and PRS's 
cog <- left_join(cog,cog_PRS, by = "Sample_Name")
cog <- left_join(cog,AD_PRS, by = "Sample_Name")

#subset to complete cases 
cog <- cog[complete.cases(cog),]

#check dimensions
dim(cog)
#[1] 19313     7

#read in protein data
prot <- readRDS("GS_ProteinGroups_RankTransformed_23Aug2023.rds")



#read in covars and select columns of interest
smid <- read.csv("2023-02-13_simd.csv")

bmi <- read.csv("2023-02-15_bmi.csv")

agesex<- read.csv("2023-03-17_agesex.csv")

alc <- read.csv("2023-02-01_alc.csv")

covars <- read.csv("2024-11-13_phenotypes.csv")

covars <- covars[,c(1,4)]
names(covars)[1] <- "id"

covars2 <- read.csv("2024-10-30_phenotypes.csv")

covars2 <- covars2[,c(1,2,4,5)]

#load libraries 
library(dplyr)
library(purrr)

# List of all data frames to be joined
dfs <- list(smid, bmi, alc, agesex, covars, covars2)

# Perform full join on all data frames by 'id'
df <- reduce(dfs, full_join, by = "id")

#drop the usual column from alc dataset 
df <- df %>% select(-usual)

#read in targets file 
targets <- readRDS("GS20k_Targets_18869.rds")

#select id and sentrix id 
targets <- targets %>% select(Sample_Name, X)

#rename columns 
names(targets) <- c("id", "SampleName")

#read in epismoker 
epismoke <- readRDS("EpiSmokEr_18859.rds")

#join targets info to epismoker
epismoke <- left_join(epismoke, targets, by = "SampleName")

#change id to numeric 
epismoke$id <- as.numeric(epismoke$id)

#merge epismoker with other covariates 
df <- full_join(df, epismoke, by = "id")

#check columns
names(df)
# [1] "id"           "rank"         "bmi"          "units"        "age"
# [6] "sex"          "depression_Y" "years"        "diabetes_Y"   "high_BP_Y"
# [11] "SampleName"   "smokingScore"

#drop methylation sample name
df <- df %>% select(-SampleName)

#read in apoe status 
apoe <- read.csv("apoe_status.csv")
apoe <- apoe[,c(1,3)]
names(apoe)[1] <- "id"

#join apoe with rest of covars
df <- full_join(df, apoe, by = "id")


#subset covars to protein and cog data 
df2 <- df[df$id %in% cog$Sample_Name,]
df2 <- df2[df2$id %in% prot$id, ]

dim(df)
#[1] 24087    12
dim(df2)
#[1] 14537    12



#create a list of covars without BMI
variables <- c("units", "rank", "smokingScore", "years",
               "age", "sex", "depression_Y",
               "diabetes_Y", "high_BP_Y")

# Efficient NA count for each variable
sapply(df2[variables], function(x) sum(is.na(x)))
# units         rank smokingScore        years          age          sex
# 1315          974          453          815           81           81
# depression_Y   diabetes_Y    high_BP_Y
# 320          323          321


remove_outliers <- function(df, columns) {
  # Identify rows to remove
  outlier_rows <- rep(FALSE, nrow(df))  # Ensure logical vector starts as FALSE
  
  for (col in columns) {
    # Calculate the mean and standard deviation of the column
    mean_val <- mean(df[[col]], na.rm = TRUE)
    sd_val <- sd(df[[col]], na.rm = TRUE)
    
    # Identify outliers (ignoring NA values)
    outliers <- (!is.na(df[[col]])) & ((df[[col]] > (mean_val + 3.5 * sd_val)) | (df[[col]] < (mean_val - 3.5 * sd_val)))
    
    # Print number of outliers in the column
    cat("Column", col, "has", sum(outliers, na.rm = TRUE), "outliers.\n")
    
    # Combine with previous outliers
    outlier_rows <- outlier_rows | outliers
  }
  
  # Remove rows that contain outliers
  df <- df[!outlier_rows, ]
  
  return(df)
}



#apply outlier removal function to columns 
df3 <- remove_outliers(df2, c("units", "rank", "smokingScore", "years"))
# Column units has 136 outliers.
# Column rank has 0 outliers.
# Column smokingScore has 85 outliers.
# Column years has 0 outliers.

dim(df2)
#[1] 14537    12
dim(df3)
#[1] 14321    12

#compare number of NAs 
sapply(df2[variables], function(x) sum(is.na(x)))
# units         rank smokingScore        years          age          sex
# 1315          974          453          815           81           81
# depression_Y   diabetes_Y    high_BP_Y
# 320          323          321

sapply(df3[variables], function(x) sum(is.na(x)))
# units         rank smokingScore        years          age          sex
# 1306          954          444          798           81           81
# depression_Y   diabetes_Y    high_BP_Y
# 316          319          317



table(is.na(df3$bmi))

# FALSE  TRUE
# 14240    80


#remove bmi less than 17 and greater than 50
df3 <- df3[is.na(df3$bmi) | (df3$bmi >= 17 & df3$bmi <= 50), ]

table(is.na(df3$bmi))

# FALSE  TRUE
# 14178    80

#create a list of covars with BMI
variables_BMI <- c("units", "rank", "smokingScore", "years",
               "age", "sex", "depression_Y",
               "diabetes_Y", "high_BP_Y", "bmi")

sapply(df3[variables_BMI], function(x) sum(is.na(x)))
# units         rank smokingScore        years          age          sex
# 1242          872          362          725            0            0
# depression_Y   diabetes_Y    high_BP_Y          bmi
# 239          242          240           80

dim(df3)
#[1] 14258    12


sapply(df3[variables_BMI], function(x) (sum(is.na(x))/14258)*100)
# units         rank smokingScore        years          age          sex
# 8.7108991    6.1158648    2.5389255    5.0848646    0.0000000    0.0000000
# depression_Y   diabetes_Y    high_BP_Y          bmi
# 1.6762519    1.6972927    1.6832655    0.5610885


#set seed 
set.seed(1234)

#convert variables to factors 
cols_to_factor <- c("sex", "depression_Y", "diabetes_Y", "high_BP_Y", "e4_count")
df3[cols_to_factor] <- lapply(df3[cols_to_factor], as.factor)

#max 60% missing 
max_missing <- 0.6

#check with rows with more than 60% missing 
rows_with_high_missing <- which(rowMeans(is.na(df3)) > max_missing)

length(rows_with_high_missing)
#[1] 0


ids <- df3$id

#remove id column
df3 <- df3[,-1]


#scale numeric 
df3_scaled <- df3  # Copy original df3
df3_scaled[, sapply(df3, is.numeric)] <- scale(df3[, sapply(df3, is.numeric)])


# rank              bmi            units     #load VIM
library(VIM)

#KNN impute missing covars
covf3 <- kNN(df3_scaled, k = 10)


sapply(covf3, function(y) sum(length(which(is.na(y)))))
      

covf3$id <- ids

dim(covf3)
#[1] 14258    23

#write out files 
setwd( "")
write.csv(covf3, file = "covars_for_cog_post_outlier_removal_KNNimputed.csv", row.names = F)
