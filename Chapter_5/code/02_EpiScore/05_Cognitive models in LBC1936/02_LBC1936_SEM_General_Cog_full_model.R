#### library ####
library("dplyr")
library("lavaan")
library("foreach")
library("doParallel")
library("readr")

## read in LBC1936 data
Lipid <- read.csv("LBC1936_with_all_lipid_EpiScores_array_date_adj.csv")
cog <- read.csv("LBC1936_WAVE1_cohort_data_KORA_EpiScore_Adjusted.csv")

#select cognitive test
cog <- select(cog, -bmi_w1, -sex)

#select episcores, measured traits, age, sex
Lipid <- Lipid[,c(4,7, 9,11, 110:115)]

#merge datasets together
data <- cog %>% full_join(Lipid, by = "lbc36no")

#set sex as factor
data$sex <- as.factor(data$sex)

#read in epismoker 
epismoke <- readRDS("lbc_epismoker.rds") 
names(epismoke)[1] <- "Basename"

#merge in epismoker
data <- data %>% left_join(epismoke, by = "Basename")

#log transform bmi and scale other traits and episcores
data$bmi_w1 <- log(data$bmi_w1)
data$bld_choles_w1 <- as.vector(scale(data$bld_choles_w1))
data$bld_hdlchol_w1 <- as.vector(scale(data$bld_hdlchol_w1))
data$bld_BMI_Epi <- as.vector(scale(data$BMI_Epi))
data$bld_HDL_Epi <- as.vector(scale(data$HDL_Epi))
data$bld_Body_fat_Epi <- as.vector(scale(data$Body_fat_Epi))
data$bld_Total_chol_Epi <- as.vector(scale(data$Total_chol_Epi))
data$bld_WHR_Epi <- as.vector(scale(data$WHR_Epi))
data$bld_Glucose_Epi <- as.vector(scale(data$Glucose_Epi))

#scale cog tests to aid in model convergance
dset_mod <- mutate(data,
                   blkdes_w1 = blkdes_w1/2,
                   blkdes_w2 = blkdes_w2/2,
                   blkdes_w3 = blkdes_w3/2,
                   blkdes_w4 = blkdes_w4/2,
                   blkdes_w5 = blkdes_w5/2,
                   vftot_w1 = vftot_w1/2,
                   vftot_w2 = vftot_w2/2,
                   vftot_w3 = vftot_w3/2,
                   vftot_w4 = vftot_w4/2,
                   vftot_w5 = vftot_w5/2,
                   lmtotal_w1 = lmtotal_w1/3,
                   lmtotal_w2 = lmtotal_w2/3,
                   lmtotal_w3 = lmtotal_w3/3,
                   lmtotal_w4 = lmtotal_w4/3,
                   lmtotal_w5 = lmtotal_w5/3,
                   digback_w1 = 3*digback_w1,
                   digback_w2 = 3*digback_w2,
                   digback_w3 = 3*digback_w3,
                   digback_w4 = 3*digback_w4,
                   digback_w5 = 3*digback_w5,
                   digsym_w1 = digsym_w1/2,
                   digsym_w2 = digsym_w2/2,
                   digsym_w3 = digsym_w3/2,
                   digsym_w4 = digsym_w4/2,
                   digsym_w5 = digsym_w5/2,
                   ittotal_w1 = ittotal_w1/2,
                   ittotal_w2 = ittotal_w2/2,
                   ittotal_w3 = ittotal_w3/2,
                   ittotal_w4 = ittotal_w4/2,
                   ittotal_w5 = ittotal_w5/2,
                   crtmean_w1 = -50 * crtmean_w1,
                   crtmean_w2 = -50 * crtmean_w2,
                   crtmean_w3 = -50 * crtmean_w3,
                   crtmean_w4 = -50 * crtmean_w4,
                   crtmean_w5 = -50 * crtmean_w5)

#summarise deprivation
summary(dset_mod$depind_w1)
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
#4    3134    5354    5265    6264   99999

#recoding spss NA code (99999) to "NA" 
dset_mod$depind_w1 <- ifelse(dset_mod$depind_w1 == 99999, NA, dset_mod$depind_w1)

#check it worked
summary(dset_mod$depind_w1)
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's
# 4    3092    5342    4565    6258    6505       8

#check na's
table(is.na(dset_mod$bmi_w1))

#FALSE  TRUE
#860   231


table(is.na(dset_mod$smokingScore))

#FALSE  TRUE
#894   197


table(is.na(dset_mod$depind_w1))
#FALSE  TRUE
#1083     8


table(is.na(dset_mod$alcunitwk_w1))

#FALSE
#1091

table(is.na(dset_mod$yrsedu_w1))

#FALSE
#1091


#scale covariates 
dset_mod$depind_w1 <- scale(dset_mod$depind_w1)
dset_mod$smokingScore <- scale(dset_mod$smokingScore)
dset_mod$alcunitwk_w1 <- scale(dset_mod$alcunitwk_w1)
dset_mod$Age_w1 <- scale(dset_mod$Age_w1)

dset_mod$date <- as.factor(dset_mod$date)
dset_mod$array <- as.factor(dset_mod$array)



#### growth curve models for individual cognitive tests ####

pgmodel1 <- '
Imatreas =~ 1*matreas_w1 + 1*matreas_w2 + 1*matreas_w3 + 1*matreas_w4 + 1*matreas_w5
Smatreas =~ 0*matreas_w1 + 2.98*matreas_w2 + 6.75*matreas_w3 + 9.82*matreas_w4 + 12.54*matreas_w5
'
fit1 <- growth(pgmodel1, dset_mod, missing = "ml.x")
summary(fit1, standardized = T)

pgmodel2 <-'
Iblkdes =~ 1*blkdes_w1 + 1*blkdes_w2 + 1*blkdes_w3 + 1*blkdes_w4 + 1*blkdes_w5
Sblkdes=~ 0*blkdes_w1 + 2.98*blkdes_w2 + 6.75*blkdes_w3 + 9.82*blkdes_w4 + 12.54*blkdes_w5
'
fit2 <- growth(pgmodel2, dset_mod, missing = "ml.x")
summary(fit2, standardized = T)

pgmodel3 <- '
Ispantot =~ 1*spantot_w1 + 1*spantot_w2 + 1*spantot_w3 + 1*spantot_w4 + 1*spantot_w5
Sspantot=~ 0*spantot_w1 + 2.98*spantot_w2 + 6.75*spantot_w3 + 9.82*spantot_w4 + 12.54*spantot_w5
'
fit3 <- growth(pgmodel3, dset_mod, missing = "ml.x")
summary(fit3, standardized = T)

pgmodel4 <- '
Inart =~ 1*nart_w1 + 1*nart_w2 + 1*nart_total_w3 + 1*nart_total_w4 + 1*nart_total_w5
Snart =~ 0*nart_w1 + 2.98*nart_w2 + 6.75*nart_total_w3 + 9.82*nart_total_w4 + 12.54*nart_total_w5
'
fit4 <- growth(pgmodel4, dset_mod, missing = "ml.x")
summary(fit4, standardized = T)

pgmodel5 <- '
Iwtar =~ 1*wtar_w1 + 1*wtar_w2 + 1*wtar_total_w3 + 1*wtar_total_w4 + 1*wtar_total_w5
Swtar =~ 0*wtar_w1 + 2.98*wtar_w2 + 6.75*wtar_total_w3 + 9.82*wtar_total_w4 + 12.54*wtar_total_w5
'
fit5 <- growth(pgmodel5, dset_mod, missing = "ml.x")
summary(fit5, standardized = T)

pgmodel6 <- '
Ivftot =~ 1*vftot_w1 + 1*vftot_w2 + 1*vftot_w3 + 1*vftot_w4 + 1*vftot_w5
Svftot =~ 0*vftot_w1 + 2.98*vftot_w2 + 6.75*vftot_w3 + 9.82*vftot_w4 + 12.54*vftot_w5
'
fit6 <- growth(pgmodel6, dset_mod, missing = "ml.x")
summary(fit6, standardized = T)

pgmodel7 <-'
Ivpatotal =~ 1*vpatotal_w1 + 1*vpatotal_w2 + 1*vpatotal_w3 + 1*vpatotal_w4 + 1*vpa_total_w5
Svpatotal =~ 0*vpatotal_w1 + 2.98*vpatotal_w2 + 6.75*vpatotal_w3 + 9.82*vpatotal_w4 + 12.54*vpa_total_w5
'
fit7 <- growth(pgmodel7, dset_mod, missing = "ml.x")
summary(fit7, standardized = T)

pgmodel8 <- '
Ilmtotal =~ 1*lmtotal_w1 + 1*lmtotal_w2 + 1*lmtotal_w3 + 1*lmtotal_w4 + 1*lmtotal_w5
Slmtotal =~ 0*lmtotal_w1 + 2.98*lmtotal_w2 + 6.75*lmtotal_w3 + 9.82*lmtotal_w4 + 12.54*lmtotal_w5
'
fit8 <- growth(pgmodel8, dset_mod, missing = "ml.x")
summary(fit8, standardized = T)

pgmodel9 <- '
Idigback =~ 1*digback_w1 + 1*digback_w2 + 1*digback_w3 + 1*digback_w4 + 1*digback_w5
Sdigback =~ 0*digback_w1 + 2.98*digback_w2 + 6.75*digback_w3 + 9.82*digback_w4 + 12.54*digback_w5
'
fit9 <- growth(pgmodel9, dset_mod, missing = "ml.x")
summary(fit9, standardized = T)

pgmodel10 <- '
Isymsear =~ 1*symsear_w1 + 1*symsear_w2 + 1*symsear_w3 + 1*symsear_w4 + 1*symsear_w5
Ssymsear =~ 0*symsear_w1 + 2.98*symsear_w2 + 6.75*symsear_w3 + 9.82*symsear_w4 + 12.54*symsear_w5
'
fit10 <- growth(pgmodel10, dset_mod, missing = "ml.x")
summary(fit10, standardized = T)

pgmodel11 <- '
Idigsym =~ 1*digsym_w1 + 1*digsym_w2 + 1*digsym_w3 + 1*digsym_w4 + 1*digsym_w5
Sdigsym =~ 0*digsym_w1 + 2.98*digsym_w2 + 6.75*digsym_w3 + 9.82*digsym_w4 + 12.54*digsym_w5
'
fit11 <- growth(pgmodel11, dset_mod, missing = "ml.x")
summary(fit11, standardized = T)

pgmodel12 <- '
Iittotal =~ 1*ittotal_w1 + 1*ittotal_w2 + 1*ittotal_w3 + 1*ittotal_w4 + 1*ittotal_w5
Sittotal =~ 0*ittotal_w1 + 2.98*ittotal_w2 + 6.75*ittotal_w3 + 9.82*ittotal_w4 + 12.54*ittotal_w5
'
fit12 <- growth(pgmodel12, dset_mod, missing = "ml.x")
summary(fit12, standardized = T)

pgmodel13 <- '
Icrtmean =~ 1*crtmean_w1 + 1*crtmean_w2 + 1*crtmean_w3 + 1*crtmean_w4 + 1*crtmean_w5
Scrtmean =~ 0*crtmean_w1 + 2.98*crtmean_w2 + 6.75*crtmean_w3 + 9.82*crtmean_w4 + 12.54*crtmean_w5
'
fit13 <- growth(pgmodel13, dset_mod, missing = "ml.x")
summary(fit13, standardized = T)



general_4p <- 'Ig =~  Iblkdes + Imatreas  + Ispantot + Ivftot + Ivpatotal + Ilmtotal +
  Idigback + Isymsear + Idigsym + Iittotal + Icrtmean + Inart + Iwtar

Sg =~ Sblkdes + Smatreas + Sspantot + Svftot + Svpatotal + Slmtotal +
  Sdigback + Ssymsear + Sdigsym + Sittotal + Scrtmean + Snart + Swtar

#indicator as scaling reference: loading=1, int=0
Iblkdes ~ 0*1
Sblkdes ~ 0*1 

# within-wave covariances between nart and wtar
nart_w1 ~~ wtar_w1
nart_w2 ~~ wtar_w2
nart_total_w3 ~~ wtar_total_w3
nart_total_w4 ~~ wtar_total_w4
nart_total_w5 ~~ wtar_total_w5

# within-test intercept-slope covariances
Imatreas ~~ Smatreas
Iblkdes ~~ Sblkdes
#Ispantot ~~Sspantot
Inart ~~ Snart
Iwtar ~~ Swtar
Ivftot ~~ Svftot
Ivpatotal ~~ Svpatotal
Ilmtotal ~~ Slmtotal
Idigback ~~ Sdigback
Isymsear ~~ Ssymsear
Idigsym ~~ Sdigsym
Iittotal ~~ Sittotal
Icrtmean ~~ Scrtmean


# within-domain intercept-intercept and slope-slope covariances
Iblkdes ~~ Imatreas # Visuospatial domain
Iblkdes ~~ Ispantot
Imatreas ~~ Ispantot
Sblkdes ~~ Smatreas 
#Sblkdes ~~ Sspantot
#Smatreas ~~ Sspantot

Inart ~~ Ivftot #Crystalized domain
Iwtar ~~ Ivftot
Iwtar ~~ Inart

Snart ~~ Svftot
Swtar ~~ Svftot
Swtar ~~ Snart

Ilmtotal ~~ Ivpatotal # Verbal memory domain
Ilmtotal ~~ Idigback
Ivpatotal ~~ Idigback
Slmtotal ~~ Svpatotal
Slmtotal ~~ Sdigback
Svpatotal ~~ Sdigback

Iittotal ~~ Idigsym #Processing speed domain
Iittotal ~~ Isymsear
Iittotal ~~ Icrtmean
Idigsym ~~ Isymsear
Idigsym ~~ Icrtmean
Isymsear ~~ Icrtmean
Sittotal ~~ Sdigsym 
Sittotal ~~ Ssymsear
Sittotal ~~ Scrtmean
Sdigsym ~~ Ssymsear
Sdigsym ~~ Scrtmean
Ssymsear ~~ Scrtmean

#fixed negative residual variance to 0 
Sspantot ~~ 0*Sspantot
'


#perform growth curve model
fitGen_4p <- growth(model = c(pgmodel1, pgmodel2, pgmodel3, pgmodel4, pgmodel5, pgmodel6, pgmodel7, pgmodel8, pgmodel9, pgmodel10, pgmodel11, pgmodel12, pgmodel13, general_4p), dset_mod,  missing = "ml.x")

#extract fit measures
fitmeasures(fitGen_4p, c("cfi", "tli", "RMSEA", "SRMR"))
#cfi   tli rmsea  srmr
#0.960 0.959 0.028 0.061

#summarise growth model 
summary(fitGen_4p, standardized = T)


################################################################
#identify and registers core for running analysis in paralell
cores <- detectCores()
cl <- makeCluster(9)
registerDoParallel(cl)

#make list of episcores 
list_e <- colnames(dset_mod[434:442])
episcores = colnames(dset_mod)[which(colnames(dset_mod)%in% list_e)]

#loop that regresses episcores or measured traits (tmp) on general cognitive function level (Ig)
foreach(i = episcores, .packages = c("lavaan")) %dopar% { 
  
  dset_mod$tmp = dset_mod[,i]
  
  reg_Ig <- '
 Ig ~ tmp + Age_w1 + sex + smokingScore + depind_w1 + alcunitwk_w1 
'
  
  fitreg_g <- growth(model = c(pgmodel1, pgmodel2, pgmodel3, pgmodel4, pgmodel5, pgmodel6, pgmodel7, pgmodel8, pgmodel9, pgmodel10, pgmodel11, pgmodel12, pgmodel13, general_4p, reg_Ig), dset_mod,  missing = "ml.x")
  
  #extract beta, se, p-value, CI's and fit measures
  output = standardizedSolution(fitreg_g)
  ind = which(output$lhs == "Ig" & output$op == "~" & output$rhs=="tmp")
  Beta = output[ind, "est.std"] 
  SE <- output[ind, "se"]
  p <- output[ind, "pvalue"]
  n <- nobs(fitreg_g)
  ci.upper <- output[ind, "ci.upper"]
  ci.lower <- output[ind, "ci.lower"]
  fitmeasures <- fitMeasures(fitreg_g, c("cfi","rmsea","srmr", "tli"))
  cfi <- fitmeasures["cfi"]
  rmsea <- fitmeasures["rmsea"]
  srmr <- fitmeasures["srmr"]
  tli <- fitmeasures["tli"]
  
  #make empty dataframe to write results to
  results <- data.frame(SeqId = NA, n = NA, beta = NA , SE = NA, P = NA, ci.upper = NA, ci.lower = NA, cfi = NA, rmsea = NA, srmr = NA, tli = NA)
  
  
  results[1,1] <- i
  results[1,2] <- n
  results[1,3] <- Beta
  results[1,4] <- SE
  results[1,5] <- p
  results[1,6] <- ci.upper
  results[1,7] <- ci.lower
  results[1,8] <- cfi
  results[1,9] <- rmsea
  results[1,10] <- srmr
  results[1,11] <- tli
  
  #save out results
  write.csv(results, paste0(i, "_G_intercept_No_Domain_LBC1936_assocs_lipid_EpiScores_fullmodel_no_yrsedu.csv"), row.names = F)
  
  
}

#loop that regresses episcores or measured traits (tmp) on general cognitive function change (Sg)
foreach(i = episcores, .packages = c("lavaan")) %dopar% { 
  
  dset_mod$tmp = dset_mod[,i]
  
  reg_Sg <- '
 Sg ~ tmp + Age_w1 + sex + smokingScore + depind_w1 + alcunitwk_w1 
'
  
  fitreg_g <- growth(model = c(pgmodel1, pgmodel2, pgmodel3, pgmodel4, pgmodel5, pgmodel6, pgmodel7, pgmodel8, pgmodel9, pgmodel10, pgmodel11, pgmodel12, pgmodel13, general_4p, reg_Sg), dset_mod,  missing = "ml.x")
  
  #extract beta, se, p-value, CI's and fit measures 
  output = standardizedSolution(fitreg_g)
  ind = which(output$lhs == "Sg" & output$op == "~" & output$rhs=="tmp")
  Beta = output[ind, "est.std"] 
  SE <- output[ind, "se"]
  p <- output[ind, "pvalue"]
  n <- nobs(fitreg_g)
  ci.upper <- output[ind, "ci.upper"]
  ci.lower <- output[ind, "ci.lower"]
  fitmeasures <- fitMeasures(fitreg_g, c("cfi","rmsea","srmr", "tli"))
  cfi <- fitmeasures["cfi"]
  rmsea <- fitmeasures["rmsea"]
  srmr <- fitmeasures["srmr"]
  tli <- fitmeasures["tli"]
  
  #create empty dataframe to write results too 
  results <- data.frame(SeqId = NA, n = NA, beta = NA , SE = NA, P = NA, ci.upper = NA, ci.lower = NA, cfi = NA, rmsea = NA, srmr = NA, tli = NA)
  
  
  results[1,1] <- i
  results[1,2] <- n
  results[1,3] <- Beta
  results[1,4] <- SE
  results[1,5] <- p
  results[1,6] <- ci.upper
  results[1,7] <- ci.lower
  results[1,8] <- cfi
  results[1,9] <- rmsea
  results[1,10] <- srmr
  results[1,11] <- tli
  
  #save out results
  write.csv(results, paste0(i, "_G_slope_No_Domain_LBC1936_assocs_lipid_EpiScores_fullmodel_no_yrsedu.csv"), row.names = F)
  
  
}

