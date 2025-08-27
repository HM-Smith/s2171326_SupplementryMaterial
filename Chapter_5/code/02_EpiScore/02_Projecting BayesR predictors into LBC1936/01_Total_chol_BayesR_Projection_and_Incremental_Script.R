###########################################################################################
######## Script to generate DNAm predictors of traits in LBC1936 at wave 1 (age 70)  ######
########################################################################################### 

### read in weights ###   
g_wts <- read.table("Bayes_Results/Total_cholesterol_prediction.txt", header = T)
meanBetas <- g_wts[,1:2]

### subset DNAm object to LBC36 wave 1 data ###
tmp36 <- which(rownames(meth1) %in% d36_w1$Basename)
meth36 <- meth1[tmp36,]

### subset to relevant CpG sites ###
a36 = which(names(meth36) %in% meanBetas$CpG)
meth36a <- meth36[,a36]

### replace missing values with mean ###
library(missMethods)
meth36b <- impute_mean(meth36a, type = "columnwise") #dim(meth36b) - [1]  861 1905

##convert dataframe to matrix ## 
meth36b <- as.matrix(meth36b) # dim(meth36b) - [1]  861 1905



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

### correlation and variance explained ###
r <- cor(ph36$bld_choles_w1, ph36$pred, use="pairwise.complete.obs")

r
#[1] 0.2312624





### Incremental R2
null <- summary(lm(bld_choles_w1 ~ age + sex, data=ph36))$r.squared
full <- summary(lm(bld_choles_w1 ~ age + sex + pred, data=ph36))$r.squared


null * 100
#[1] 11.4182

full * 100
#[1] 14.97286



round(100*(full - null), 3)
#[1] 3.555




#drop levels that are redundant 
ph36$date <- droplevels(ph36$date)


### Incremental R2
null <- summary(lm(bld_choles_w1 ~ age + sex + as.factor(date) + as.factor(array), data=ph36))$r.squared
full <- summary(lm(bld_choles_w1 ~ age + sex + as.factor(date) + as.factor(array) + pred, data=ph36))$r.squared


round(100*(full - null), 3)
#[1] 3.194



names(ph36)[107] <- "Total_chol_Epi"

write.csv(ph36, file = "LBC1936_with_Total_cholesterol_BayesR.csv")