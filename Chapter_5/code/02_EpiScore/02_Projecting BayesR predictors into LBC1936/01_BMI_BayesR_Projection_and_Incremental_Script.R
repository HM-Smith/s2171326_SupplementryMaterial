###########################################################################################
######## Script to generate DNAm predictors of traits in LBC1936 at wave 1 (age 70)  ######
########################################################################################### 

# Load requisite libraries 
library(data.table)


### load LBC target file ##
d = readRDS("targets_3489_bloodonly.rds")

#subset to  LBC1936 wave 1 and set 1 
d36_w1 <- d[d$cohort=="LBC36" & d$WAVE==1 & d$set==1,]

### load in methylation data ###
dat <- readRDS("LBC_betas_3489_bloodonly.rds")
meth = t(dat)
meth1 = as.data.frame(meth)
meth1$id = as.character(rownames(meth1))


### read in weights ###   
g_wts <- read.table("/Bayes_Results/bmi_prediction.txt", header = T) #dim(g_wts) [1] 386399      5
dim(g_wts)
#[1] 386399      5


meanBetas <- g_wts[,1:2]

### subset DNAm object to LBC36 wave 1 data ###
tmp36 <- which(rownames(meth1) %in% d36_w1$Basename) # length(tmp36) - [1] 861

meth36 <- meth1[tmp36,] # dim(meth36) - [1]    861 459310


### subset to relevant CpG sites ###
a36 = which(names(meth36) %in% meanBetas$CpG) #length(a36) - [1] 374995

meth36a <- meth36[,a36] #dim(meth36a) - [1]  861 374995



### replace missing values with mean ###
library(missMethods)
meth36b <- impute_mean(meth36a, type = "columnwise") #dim(meth36b) - [1]    861 374995


##convert dataframe to matrix ## 
meth36b <- as.matrix(meth36b) # dim(meth36b) - [1]    861 374995


### line up weights and CpGs ###
b36 = which(meanBetas$CpG %in% colnames(meth36b))
length(b36)
#[1] 374995 - overlap between LBC1936 meth and CpGs in predictor file 
dim(meanBetas)
#[1] 386399      2  - number of CpGs in predictor file 


mean_betas <- meanBetas[b36,] #dim(mean_betas) - [1] 374995      2


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
r <- cor(log(ph36$bmi_w1), ph36$pred, use="pairwise.complete.obs")

r

#[1] 0.3577143


### Incremental R2
null <- summary(lm(log(bmi_w1) ~ age + as.factor(sex), data=ph36))$r.squared
full <- summary(lm(log(bmi_w1) ~ age + as.factor(sex) + pred, data=ph36))$r.squared

null * 100
#[1] 1.367787 

full * 100
# [1] 14.0606

round(100*(full - null), 3)
#[1] 12.693


#drop levels that are redundant 
ph36$date <- droplevels(ph36$date)


### Incremental R2 with covariates
null <- summary(lm(log(bmi_w1) ~ age + as.factor(sex) + as.factor(date) + as.factor(array), data=ph36))$r.squared
full <- summary(lm(log(bmi_w1) ~ age + as.factor(sex) + as.factor(date) + as.factor(array) + pred, data=ph36))$r.squared
round(100*(full - null), 3)
#[1] 12.354



## not enough levels 
names(ph36)[107] <- "BMI_Epi"
write.csv(ph36, file = "LBC1936_with_BMI_EpiScore_BayesR.csv")




