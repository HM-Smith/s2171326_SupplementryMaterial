###########################################################################################
######## Script to generate DNAm predictors of traits in LBC1936 at wave 1 (age 70)  ######
########################################################################################### 

########################## ELNET PREDICTOR OF HDL #########################################

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
g_wts <- readRDS("HDL_cholesterol_predictor.rds")
dim(g_wts)
#[1] 1676    3

meanBetas <- g_wts[,2:1]
names(meanBetas)=c("CpG","Beta")

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
length(b36)
#[1] 1633
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
r <- cor(ph36$bld_hdlchol_w1, ph36$pred, use="pairwise.complete.obs")

r
#[1] 0.4808707



### Incremental R2
null <- summary(lm(bld_hdlchol_w1 ~ age + sex, data=ph36))$r.squared
full <- summary(lm(bld_hdlchol_w1 ~ age + sex + pred, data=ph36))$r.squared

null * 100
#[1] 9.893691

  

full * 100
#[1] 29.19232



round(100*(full - null), 3)
#[1] 19.299



ph36$date <- droplevels(ph36$date)


### Incremental R2 with covariates
null <- summary(lm(bld_hdlchol_w1 ~ age + as.factor(sex) + as.factor(date) + as.factor(array), data=ph36))$r.squared
full <- summary(lm(bld_hdlchol_w1 ~ age + as.factor(sex) + as.factor(date) + as.factor(array) + pred, data=ph36))$r.squared
round(100*(full - null), 3)
#[1] 18.515



names(ph36)[107] <- "HDL_Epi"

write.csv(ph36, file = "LBC1936_with_HDL.csv")


