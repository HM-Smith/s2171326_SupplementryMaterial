###########################################################################################
######## Script to generate DNAm predictors of traits in LBC1936 at wave 1 (age 70)  ######
########################################################################################### 

# Load requisite libraries 
library(data.table)


### load LBC target file ###   # This just does Wave 1 
d = readRDS("targets_3489_bloodonly.rds")
d36_w1 <- d[d$cohort=="LBC36" & d$WAVE==1 & d$set==1,]

### load in methylation data ###
dat <- readRDS("LBC_betas_3489_bloodonly.rds")
meth = t(dat)
meth1 = as.data.frame(meth)
meth1$id = as.character(rownames(meth1))


### read in weights ###   
g_wts <- read.table("Bayes_Results/body_fat_prediction.txt", header = T)
dim(g_wts)
#[1] 386399      5

meanBetas <- g_wts[,1:2]

### subset DNAm object to LBC36 wave 1 data ###
tmp36 <- which(rownames(meth1) %in% d36_w1$Basename)
meth36 <- meth1[tmp36,]

### subset to relevant CpG sites ###
a36 = which(names(meth36) %in% meanBetas$CpG)
meth36a <- meth36[,a36]


### replace missing values with mean ###
library(missMethods)
meth36b <- impute_mean(meth36a, type = "columnwise") #dim(meth36b) - [1]    861 374995


##convert dataframe to matrix ## 
meth36b <- as.matrix(meth36b) 

### line up weights and CpGs ###
b36 = which(meanBetas$CpG %in% colnames(meth36b)) 
mean_betas <- meanBetas[b36,]
meth36c <- meth36b[,match(colnames(meth36b), mean_betas$CpG)]
meth36d <- meth36c[,mean_betas$CpG]

### create predictor ###
pred36 <- meth36d %*% mean_betas$Beta
pred36 <- as.data.frame(pred36)
names(pred36) <- "pred"
pred36$Basename <- rownames(pred36)
d36 <- merge(d36_w1, pred36, by="Basename")


### load in phenotypes ###
library("foreign")
ph36 <- read.spss("EpigeneticLipidScoresAndBrainHealth_HS_27FEB2023.sav", to.data.frame=T, use.value.labels=F)
ph36$sex=NULL # sex is already in other dataframe so not needed twice 

### merge in phenotype data with predictor 
ph36=merge(ph36,d36,by.x="lbc36no",by="ID")


names(ph36)[107] <- "Body_fat_Epi"

write.csv(ph36, file = "LBC1936_with_Body_fat_EpiScore_BayesR.csv")