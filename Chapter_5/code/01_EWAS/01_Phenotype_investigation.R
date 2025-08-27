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



ped1 <- read.csv("pedigree.csv")
ped2 = data.frame(famid=c(,), volid=c(,), father=c(0,0), mother=c(0,0), sex= c("F", "F"))
ped = rbind(ped1, ped2)
kin <- with(ped, pedigree(volid, father, mother, sex, famid=famid))
kin_model <- kinship(kin) 



## Subset cov file to those in pedigree file to avoid error in lmekin function (can't have anyone in phenotype file thats not in kinship matrix)
cov = combined_data[which(combined_data$ID %in% ped$volid),]

#check dimensions of data 
dim(combined_data)
#[1] 18413    20


dim(cov)
#[1] 18411    20



####### BMI ###########


## filter out NA's 
#bmi_filtered <- cov %>% filter(bmi != "NA") 
bmi_filtered <- cov

#check data dimensions 
dim(bmi_filtered)
#[1] 18411    20

#dat summary of bmi
summary(bmi_filtered$bmi)
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's
#10.49   23.04   25.91   26.67   29.26   71.35     114

#filter for BMI >= 17 and <= 50
bmi_update <- bmi_filtered[bmi_filtered$bmi >= 17 & bmi_filtered$bmi <= 50, ]

# check data dimensions 
dim(bmi_update)
#[1] 18329    20


#summary data of BMI
summary(bmi_update$bmi)
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's
#17.00   23.06   25.91   26.64   29.25   49.99     114

# plot hist of BMI
hist(bmi_update$bmi, breaks = 100)

#plot BMI against body fat
plot(bmi_update$bmi, bmi_update$body_fat)

#summarise body fat
summary(bmi_update$body_fat)
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's
#1.0    23.0    29.7    29.9    37.0    62.0     494

#filter body fat >= 8 and <= 50
body_fat_update <- bmi_update[bmi_update$body_fat >= 8 & bmi_update$body_fat <= 50,]

# check dimensions of data
dim(body_fat_update)
#[1] 17993    20

#summary of body fat
summary(body_fat_update$body_fat)
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's
#8.00   23.10   29.70   29.84   36.80   50.00     494

#plot bmi against body fat in new filtered file
plot(body_fat_update$bmi, body_fat_update$body_fat)

#plot bmi against whr
plot(body_fat_update$bmi, body_fat_update$whr)

#filter data by whr < 1.5
whr_update <- body_fat_update[body_fat_update$whr < 1.5, ]

#check data dimensions 
dim(whr_update)
#[1] 17989    20

#plot bmi against whr
plot(whr_update$bmi, whr_update$whr)

#plot bmi against whr
plot(whr_update$bmi, whr_update$body_fat)

#manual check and removal of outliers based on plots
check <- whr_update[whr_update$bmi < 18 & whr_update$body_fat > 35,]

check2 <- check[!is.na(check$bmi),]

check2$ID #id for indvidual with BMI = 17.9 and body fat % = 39.5
#[1] 73868

check <- whr_update[whr_update$bmi > 30 & whr_update$body_fat < 10, ]

check2 <- check[!is.na(check$bmi),]

check2$ID
#[1]  143137

check <- whr_update[whr_update$bmi > 34 & whr_update$body_fat < 17, ]

check2 <- check[!is.na(check$bmi),]

check2$ID
#[1] 14562


check <- whr_update[whr_update$bmi > 45 & whr_update$body_fat < 31, ]
check2 <- check[!is.na(check$bmi),]

check2$ID
#135019

#ids to remove based on above checks 
id_remove <- c("73868", "14562", "143137", "135019")

#removign outlier ids 
outliers_removed <- whr_update[!whr_update$ID %in% id_remove, ]

#check data dimensions 
dim(whr_update)
#[1] 17989    20
dim(outliers_removed)
#[1] 17985    20

#plot body fat, bmi and whr against eachother to do final check for outliers 
plot(outliers_removed$bmi, outliers_removed$body_fat)

plot(outliers_removed$bmi, outliers_removed$whr)

plot(outliers_removed$body_fat, outliers_removed$whr)


#################################################################################################################################################################################################################
## filter out NA's 
bmi <- outliers_removed %>% filter(bmi != "NA") 

dim(bmi)
#[1] 17304    20


## Linear mixed-effects model  
bmi_mod = lmekin(log(bmi[,"bmi"]) ~ bmi$age + bmi$age^2 + as.factor(bmi$sex)+ (1|bmi$ID), varlist = kin_model*2, na.action = na.exclude)

## store residuals 
bmi_res = as.data.frame(bmi_mod$residuals) 

## tidy up output file 
bmi_res$ID <- bmi$Sample_Sentrix_ID
bmi_res <- bmi_res[,c(2,1)]
names(bmi_res)[2] <- "bmi"

dim(bmi_res)
#[1] 17304     2

## save out filtered residuals 
write.csv(bmi_res, file = "bmi_res.csv")

#################################################################################################################################################################################################################
whr_filtered <- outliers_removed %>% filter(whr != "NA")

dim(whr_filtered)
#[1] 17304    20


## Linear mixed-effects model  
whr_mod = lmekin(whr_filtered[,"whr"] ~ whr_filtered$age + whr_filtered$age^2 + as.factor(whr_filtered$sex)  + (1|whr_filtered$ID), varlist = kin_model*2, na.action = na.exclude)

## store residuals 
whr_res = as.data.frame(whr_mod$residuals) 

## tidy up output file 
whr_res$ID <- whr_filtered$Sample_Sentrix_ID
whr_res <- whr_res[,c(2,1)]
names(whr_res)[2] <- "whr"

dim(whr_res)
#[1] 17304     2

## save out filtered residuals 
write.csv(whr_res, file = "whr_res.csv")

#################################################################################################################################################################################################################

## filter out NA's 
body_fat_filtered <- outliers_removed %>% filter(body_fat != "NA") 

dim(body_fat_filtered)
#[1] 17304    20


## Linear mixed-effects model  
body_fat_mod = lmekin(body_fat_filtered[,"body_fat"] ~ body_fat_filtered$age + body_fat_filtered$age^2 + factor(body_fat_filtered$sex)  + (1|body_fat_filtered$ID), varlist = kin_model*2, na.action = na.exclude)

## store residuals 
body_fat_res = as.data.frame(body_fat_mod$residuals) 

## tidy up output file 
body_fat_res$ID <- body_fat_filtered$Sample_Sentrix_ID
body_fat_res <- body_fat_res[,c(2,1)]
names(body_fat_res)[2] <- "body_fat"

dim(body_fat_res)
#[1] 17304     2

## save out filtered residuals 
write.csv(body_fat_res, file = "body_fat_res.csv")

#################################################################################################################################################################################################################
## filter out NA's 
Glucose <- cov %>% filter(Glucose != "NA")

dim(Glucose)
#[1] 18081    20

## dealing with ouliers
mean_Glucose <- mean(Glucose$Glucose)
sd_Glucose <- sd(Glucose$Glucose) * 4 #calculate 4 standard deviations 
upper_Glucose <- mean_Glucose + sd_Glucose #add 4 SD onto mean
lower_Glucose <- mean_Glucose - sd_Glucose # minus 4 SD from mean 

Glucose_filtered <- Glucose %>% filter(Glucose < upper_Glucose & Glucose > lower_Glucose) # filter data for upper and lower limits 

dim(Glucose_filtered)
#[1] 17908    20

## Linear mixed-effects model  
Glucose_mod = lmekin(Glucose_filtered[,"Glucose"] ~ Glucose_filtered$age + Glucose_filtered$age^2 + factor(Glucose_filtered$sex) + (1|Glucose_filtered$ID), varlist = kin_model*2, na.action = na.exclude)

## store residuals 
Glucose_res = as.data.frame(Glucose_mod$residuals) 

## tidy up output file 
Glucose_res$ID <- Glucose_filtered$Sample_Sentrix_ID
Glucose_res <- Glucose_res[,c(2,1)]
names(Glucose_res)[2] <- "Glucose"

## save out filtered residuals 
write.csv(Glucose_res, file = "Glucose_res.csv")

#################################################################################################################################################################################################################
## filter out NA's 
HDL_cholesterol <- cov %>% filter(HDL_cholesterol != "NA") 

dim(HDL_cholesterol)
#[1] 18251    20

## dealing with ouliers
mean_HDL_cholesterol <- mean(HDL_cholesterol$HDL_cholesterol)
sd_HDL_cholesterol <- sd(HDL_cholesterol$HDL_cholesterol) * 4 #calculate 4 standard deviations 
upper_HDL_cholesterol <- mean_HDL_cholesterol + sd_HDL_cholesterol #add 4 SD onto mean
lower_HDL_cholesterol <- mean_HDL_cholesterol - sd_HDL_cholesterol # minus 4  SD from mean 

HDL_cholesterol_filtered <- HDL_cholesterol %>% filter(HDL_cholesterol < upper_HDL_cholesterol & HDL_cholesterol > lower_HDL_cholesterol) # filter data for upper and lower limits 

dim(HDL_cholesterol_filtered)
#[1] 18225    20

## Linear mixed-effects model  
HDL_cholesterol_mod = lmekin(HDL_cholesterol_filtered[,"HDL_cholesterol"] ~ HDL_cholesterol_filtered$age + HDL_cholesterol_filtered$age^2 + factor(HDL_cholesterol_filtered$sex) + (1|HDL_cholesterol_filtered$ID), varlist = kin_model*2, na.action = na.exclude)

## store residuals 
HDL_cholesterol_res = as.data.frame(HDL_cholesterol_mod$residuals) 

## tidy up output file 
HDL_cholesterol_res$ID <- HDL_cholesterol_filtered$Sample_Sentrix_ID
HDL_cholesterol_res <- HDL_cholesterol_res[,c(2,1)]
names(HDL_cholesterol_res)[2] <- "HDL_cholesterol"


## save out filtered residuals 
write.csv(HDL_cholesterol_res, file = "HDL_cholesterol_res.csv")

#################################################################################################################################################################################################################
## filter out NA's 
Total_cholesterol <- cov %>% filter(Total_cholesterol != "NA") 

dim(Total_cholesterol)
#[1] 18285    20

## dealing with ouliers
mean_Total_cholesterol <- mean(Total_cholesterol$Total_cholesterol)
sd_Total_cholesterol <- sd(Total_cholesterol$Total_cholesterol) * 4 #calculate 4 standard deviations 
upper_Total_cholesterol <- mean_Total_cholesterol + sd_Total_cholesterol #add 4 SD onto mean
lower_Total_cholesterol <- mean_Total_cholesterol - sd_Total_cholesterol # minus 4 SD from mean 

Total_cholesterol_filtered <- Total_cholesterol %>% filter(Total_cholesterol < upper_Total_cholesterol & Total_cholesterol > lower_Total_cholesterol) # filter data for upper and lower limits 


## Linear mixed-effects model  
Total_cholesterol_mod = lmekin(Total_cholesterol_filtered[,"Total_cholesterol"] ~ Total_cholesterol_filtered$age + Total_cholesterol_filtered$age^2 + factor(Total_cholesterol_filtered$sex) + (1|Total_cholesterol_filtered$ID), varlist = kin_model*2, na.action = na.exclude)

## store residuals 
Total_cholesterol_res = as.data.frame(Total_cholesterol_mod$residuals) 

## tidy up output file 
Total_cholesterol_res$ID <- Total_cholesterol_filtered$Sample_Sentrix_ID
Total_cholesterol_res <- Total_cholesterol_res[,c(2,1)]
names(Total_cholesterol_res)[2] <- "Total_cholesterol"

dim(Total_cholesterol_res)
#[1] 18270     2

## save out filtered residuals 
write.csv(Total_cholesterol_res, file = "Total_cholesterol.csv")


################## save out phenotype files ##############################

pheno_bmi <- data.frame(FID = bmi_res$ID,
                        IID = bmi_res$ID,
                        bmi = scale(bmi_res$bmi))
write.table(pheno_bmi, file="bmi_res.phen", row.names=F, sep=' ')



pheno_body <- data.frame(FID = body_fat_res$ID,
                         IID = body_fat_res$ID,
                         body_fat = scale(body_fat_res$body_fat))
write.table(pheno_body, file= "body_fat_res.phen", row.names=F, sep=' ')



pheno_whr <- data.frame(FID = whr_res$ID,
                        IID = whr_res$ID,
                        whr = scale(whr_res$whr))
write.table(pheno_whr, file="whr_res.phen", row.names=F, sep=' ')


pheno_glucose <- data.frame(FID = Glucose_res$ID,
                            IID = Glucose_res$ID,
                            Glucose = scale(Glucose_res$Glucose))
write.table(pheno_glucose, file="Glucose_res.phen", row.names=F, sep=' ')



pheno_HDL <- data.frame(FID = HDL_cholesterol_res$ID,
                        IID = HDL_cholesterol_res$ID,
                        HDL = scale(HDL_cholesterol_res$HDL_cholesterol))
write.table(pheno_HDL, file="HDL_cholesterol_res.phen", row.names=F, sep=' ')



pheno_Total_chol <- data.frame(FID = Total_cholesterol_res$ID,
                               IID = Total_cholesterol_res$ID,
                               Total_cholesterol = scale(Total_cholesterol_res$Total_cholesterol))
write.table(pheno_Total_chol, file="Total_cholesterol_res.phen", row.names=F, sep=' ')


