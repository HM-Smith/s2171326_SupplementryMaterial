#load dplyr 
library(dplyr)
#read in dementia data 
time_to_event <- read.csv("time_to_event_data_dementia_death_25102024.csv")

#load data table library 
library(data.table)

#read in protein data
prot <- readRDS("GS_ProteinGroups_RankTransformed_23Aug2023.rds")


#read in covars
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
df2 <- df[df$id %in% time_to_event$Sample_Name,]
df2 <- df2[df2$id %in% prot$id, ]

dim(df)
#[1] 24087    12
dim(df2)
#[1] 7128   12


#create a list of covars without BMI
variables <- c("units", "rank", "smokingScore", "years",
               "age", "sex", "depression_Y",
               "diabetes_Y", "high_BP_Y")

# Efficient NA count for each variable
sapply(df2[variables], function(x) sum(is.na(x)))
# units         rank smokingScore        years          age          sex
# 706          440          462          357            0            0
# depression_Y   diabetes_Y    high_BP_Y
# 95           97           94



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
# Column units has 58 outliers.
# Column rank has 0 outliers.
# Column smokingScore has 37 outliers.
# Column years has 28 outliers.
  

dim(df2)
#[1] 7128   12
dim(df3)
#[1] 7008   12


#compare number of NAs 
sapply(df2[variables], function(x) sum(is.na(x)))
# units         rank smokingScore        years          age          sex
# 706          440          462          357            0            0
# depression_Y   diabetes_Y    high_BP_Y
# 95           97           94

sapply(df3[variables], function(x) sum(is.na(x)))
# units         rank smokingScore        years          age          sex
# 702          434          452          350            0            0
# depression_Y   diabetes_Y    high_BP_Y
# 94           96           93



table(is.na(df3$bmi))
# 
# FALSE  TRUE
# 6956    52



#remove BMI less than 17 and over 50
df3 <- df3[is.na(df3$bmi) | (df3$bmi >= 17 & df3$bmi <= 50), ]


dim(df3)
#[1] 6981   12

table(is.na(df3$bmi))
# 
# FALSE  TRUE
# 6929    52

#set max missing to 60%
max_missing <- 0.6

#check for missing more than 60%
rows_with_high_missing <- which(rowMeans(is.na(df3)) > max_missing)

rows_with_high_missing
# 14243 21764 24085
# 4445  6716  6981

#list of rows with high missing 
rows_with_high_missing <- names(rows_with_high_missing)


dim(df3)
#[1] 6981   12

#remove rows with high missing 
df3 <- df3[!(rownames(df3) %in% rows_with_high_missing), ]


dim(df3)
#[1] 6978   12

#create a list of covars without BMI
variables_bmi <- c("units", "rank", "smokingScore", "years",
               "age", "sex", "depression_Y",
               "diabetes_Y", "high_BP_Y", "bmi")

sapply(df3[variables_bmi], function(x) (sum(is.na(x))/6978)*100)
# units         rank smokingScore        years          age          sex
# 9.9598739    6.1478934    6.3915162    4.9297793    0.0000000    0.0000000
# depression_Y   diabetes_Y    high_BP_Y          bmi
# 1.2897678    1.3184293    1.2754371    0.7165377


#set seed 
set.seed(1234)

#convert variables to factors 
cols_to_factor <- c("sex", "depression_Y", "diabetes_Y", "high_BP_Y", "e4_count")
df3[cols_to_factor] <- lapply(df3[cols_to_factor], as.factor)

ids <- df3$id

#remove id column
df3 <- df3[,-1]


#scale numeric 
df3_scaled <- df3  # Copy original df3
df3_scaled[, sapply(df3, is.numeric)] <- scale(df3[, sapply(df3, is.numeric)])

sapply(df3[variables], function(x) sum(is.na(x)))
# units         rank smokingScore        years          age          sex
# 695          429          446          344            0            0
# depression_Y   diabetes_Y    high_BP_Y
# 90           92           89


#load VIM
library(VIM)

#KNN impute missing covars
covf3 <- kNN(df3_scaled, k = 10)

sapply(covf3, function(y) sum(length(which(is.na(y)))))
# rank              bmi            units              age
# 0                0                0                0
# sex     depression_Y            years       diabetes_Y
# 0                0                0                0
# high_BP_Y     smokingScore         e4_count         rank_imp
# 0                0                0                0
# bmi_imp        units_imp          age_imp          sex_imp
# 0                0                0                0
# depression_Y_imp        years_imp   diabetes_Y_imp    high_BP_Y_imp
# 0                0                0                0
# smokingScore_imp     e4_count_imp
# 0                0


# 
# library(ggplot2)
# 
# ggplot(covf3, aes(age, bmi, color = bmi_imp)) +
#   geom_point()

covf3$id <- ids

dim(covf3)
#[1] 6978   23


setwd( "")
write.csv(covf3, file = "covars_for_dementia_post_outlier_removal_KNNimputed.csv", row.names = F)