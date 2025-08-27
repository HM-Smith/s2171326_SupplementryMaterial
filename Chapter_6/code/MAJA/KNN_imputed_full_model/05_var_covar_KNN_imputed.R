#Load dplyr package
library(dplyr)

#read in variances/covariance 
df <- read.table("mean_V.txt")

#read in the variance of variance and covariance 
df2 <- read.table("var_V.txt")

#reintroduce column names 
names(df) <- c( "digit symbol", "verbal total", "vocabulary","logical memory", "AD PRS", "cognitive PRS")

names(df2) <- c( "digit symbol", "verbal total", "vocabulary","logical memory", "AD PRS", "cognitive PRS")

#save column names for ids 
ids <- names(df)

#extract the diagonal (variances) from df 
var <- diag(as.matrix(df))

#extract the diagonal (variance of the variance)
var_var <- diag(as.matrix(df2))

#make df with variance and variance of variance 
var_df <- data.frame(Measure = ids, Variance = var)
var_var_df <- data.frame(Measure = ids, Variance_of_variance = var_var)

#merge together 
all <- full_join(var_df, var_var_df, by = "Measure")

#multiply variance and variance of variance by 100
all$Variance <- all$Variance * 100
all$Variance_of_variance <- all$Variance_of_variance * 100

#save out variance file 
write.csv(all, file = "variance_explained.csv")

#create upper and lower variance limit 
all$lower <- all$Variance - all$Variance_of_variance
all$upper <- all$Variance + all$Variance_of_variance



#create covar df by extracting the upper triangle of dataset (lower triangle also works)
covar <- as.matrix(df)
rownames(covar) <- colnames(covar)
covariances <- covar[upper.tri(covar)]

#make row names ans column names ids 
row_names <- rownames(covar)
col_names <- colnames(covar)

#create labels for covariances 
labels <- outer(row_names, col_names, FUN = paste, sep = "-")[upper.tri(covar)]

#create dataframe with covariances and labels 
labeled_covariances <- data.frame(Covariances = covariances, Labels = labels)

#variance of covariance matrix 
covar_var <- as.matrix(df2)

#set row names and column names to be the same 
rownames(covar_var) <- colnames(covar_var)

#extract variance of covariance (upper triangle or lower triangle )
covariances_var <- covar_var[upper.tri(covar_var)]

#set row names and column names 
row_names <- rownames(covar_var)
col_names <- colnames(covar_var)

#create labels for variance of covariance matrix 
labels <- outer(row_names, col_names, FUN = paste, sep = "-")[upper.tri(covar_var)]

#make dataframe with labels and variances of covaiances 
labeled_covariances_var <- data.frame(Covariances_var = covariances_var, Labels = labels)

#merge together 
all_covar <- full_join(labeled_covariances, labeled_covariances_var, by = "Labels")

#multiple covariance and variance of covariance by 100
all_covar$Covariances <- all_covar$Covariances * 100 
all_covar$Covariances_var <- all_covar$Covariances_var * 100


#save put covariance file 
write.csv(all_covar, file = "omic_covariances.csv")
