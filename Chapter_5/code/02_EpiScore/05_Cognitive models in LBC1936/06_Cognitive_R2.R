#### library ####
library("dplyr")
library("lavaan")
library("foreach")
library("doParallel")
library("readr")

## read in LBC1936 data - episcores, measured traits, cog tests and covars
Lipid <- read.csv("LBC1936_with_all_lipid_EpiScores_array_date_adj.csv")
cog <- read.csv("LBC1936_WAVE1_cohort_data_KORA_EpiScore_Adjusted.csv")

#select columns with cog tests
cog <- select(cog, -bmi_w1, -sex)

#select episcores and measured traits and covars
Lipid <- Lipid[,c(4,7, 9,11, 110:115)]

#join datasets together
data <- cog %>% full_join(Lipid, by = "lbc36no")

#set sex as factor
data$sex <- as.factor(data$sex)

#log transform bmi
data$bmi_w1 <- log(data$bmi_w1)

#scale all variables 
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


#growth curve model for general cognitive function level and change
fitGen_4p <- growth(model = c(pgmodel1, pgmodel2, pgmodel3, pgmodel4, pgmodel5, pgmodel6, pgmodel7, pgmodel8, pgmodel9, pgmodel10, pgmodel11, pgmodel12, pgmodel13, general_4p), dset_mod,  missing = "ml.x")

#extract fit measures
fitmeasures(fitGen_4p, c("cfi", "tli", "RMSEA", "SRMR"))

#summary of growth model 
summary(fitGen_4p, standardized = T)


################################################################
#make list of episcores/measured traits
list_e <- colnames(dset_mod[435:443])

#print list
list_e
#[1] "bmi_w1"         "bld_choles_w1"  "bld_hdlchol_w1" "BMI_Epi"
#[5] "HDL_Epi"        "Total_chol_Epi" "Glucose_Epi"    "Body_fat_Epi"
#[9] "WHR_Epi"

######################################################################
## age + sex ##

#regress age and sex on cog function intercept 
reg_age_sex <- '
 Ig ~ Age_w1 + sex
'

fitreg_age_sex <- growth(model = c(pgmodel1, pgmodel2, pgmodel3, 
                               pgmodel4, pgmodel5, pgmodel6, pgmodel7, 
                               pgmodel8, pgmodel9, pgmodel10, pgmodel11, 
                               pgmodel12, pgmodel13, general_4p, reg_age_sex), 
                     dset_mod,  missing = "ml.x")

#print summary to get var explained
summary(fitreg_age_sex, standardized = T, rsquare = T)


## BMI - #regress BMI, age and sex on cog function intercept 


  reg_bmi <- '
 Ig ~ bmi_w1 + Age_w1 + sex
'
  
  fitreg_bmi <- growth(model = c(pgmodel1, pgmodel2, pgmodel3, 
                               pgmodel4, pgmodel5, pgmodel6, pgmodel7, 
                               pgmodel8, pgmodel9, pgmodel10, pgmodel11, 
                               pgmodel12, pgmodel13, general_4p, reg_bmi), 
                     dset_mod,  missing = "ml.x")
  
  #print summary to get var explained
  summary(fitreg_bmi, standardized = T, rsquare = T)

  
## BMI EpiScore - #regress BMI episcore, age and sex on cog function intercept 


  reg_Epi <- '
 Ig ~ BMI_Epi + Age_w1 + sex
'
  
  fitreg_Epi <- growth(model = c(pgmodel1, pgmodel2, pgmodel3, pgmodel4, 
                                 pgmodel5, pgmodel6, pgmodel7, pgmodel8, 
                                 pgmodel9, pgmodel10, pgmodel11, pgmodel12, 
                                 pgmodel13, general_4p, reg_Epi), 
                       dset_mod,  missing = "ml.x")
  #print summary to get var explained
  summary(fitreg_Epi, standardized = T, rsquare = T)
  

## BMI and BMI EpiScore - #regress BMI, BMI EpiScore, age and sex on cog function intercept 

  
  reg_bmi_Epi <- '
 Ig ~ bmi_w1 + BMI_Epi + Age_w1 + sex
'
  
  fitreg_bmi_Epi <- growth(model = c(pgmodel1, pgmodel2, pgmodel3, 
                               pgmodel4, pgmodel5, pgmodel6, 
                               pgmodel7, pgmodel8, pgmodel9, pgmodel10, 
                               pgmodel11, pgmodel12, pgmodel13, general_4p, reg_bmi_Epi), 
                     dset_mod,  missing = "ml.x")
  #print summary to get var explained
  summary(fitreg_bmi_Epi, standardized = T, rsquare = T)
  
#########################################################################################
  ######################################################################
  
  ## HDL - #regress HDL, age and sex on cog function intercept 

  
  reg_HDL <- '
 Ig ~ bld_hdlchol_w1 + Age_w1 + sex
'
  
  fitreg_HDL <- growth(model = c(pgmodel1, pgmodel2, pgmodel3, 
                                 pgmodel4, pgmodel5, pgmodel6, pgmodel7, 
                                 pgmodel8, pgmodel9, pgmodel10, pgmodel11, 
                                 pgmodel12, pgmodel13, general_4p, reg_HDL), 
                       dset_mod,  missing = "ml.x")
  #print summary to get var explained
  summary(fitreg_HDL, standardized = T, rsquare = T)
  
  
  ## HDL EpiScore - #regress HDL EpiScore, age and sex on cog function intercept 

  
  reg_HDLEpi <- '
 Ig ~ HDL_Epi + Age_w1 + sex
'
  
  fitreg_HDLEpi <- growth(model = c(pgmodel1, pgmodel2, pgmodel3, pgmodel4, 
                                 pgmodel5, pgmodel6, pgmodel7, pgmodel8, 
                                 pgmodel9, pgmodel10, pgmodel11, pgmodel12, 
                                 pgmodel13, general_4p, reg_HDLEpi), 
                       dset_mod,  missing = "ml.x")
  #print summary to get var explained
  summary(fitreg_HDLEpi, standardized = T, rsquare = T)
  
  
  ## HDL and HDL EpiScore - #regress HDL, HDL EpiScore, age and sex on cog function intercept 

  
  reg_HDL_HDLEpi <- '
 Ig ~ bld_hdlchol_w1 + HDL_Epi + Age_w1 + sex
'
  
  fitreg_HDL_HDLEpi <- growth(model = c(pgmodel1, pgmodel2, pgmodel3, 
                                     pgmodel4, pgmodel5, pgmodel6, 
                                     pgmodel7, pgmodel8, pgmodel9, pgmodel10, 
                                     pgmodel11, pgmodel12, pgmodel13, general_4p, 
                                     reg_HDL_HDLEpi), 
                           dset_mod,  missing = "ml.x")
  #print summary to get var explained
  summary(fitreg_HDL_HDLEpi, standardized = T, rsquare = T)
  
  
  #########################################################################################
  ######################################################################
  
  ## Total cholesterol - #regress total cholesterol, age and sex on cog function intercept 

  
  reg_total <- '
 Ig ~ bld_choles_w1 + Age_w1 + sex
'
  
  fitreg_total <- growth(model = c(pgmodel1, pgmodel2, pgmodel3, 
                                 pgmodel4, pgmodel5, pgmodel6, pgmodel7, 
                                 pgmodel8, pgmodel9, pgmodel10, pgmodel11, 
                                 pgmodel12, pgmodel13, general_4p, reg_total), 
                       dset_mod,  missing = "ml.x")
  #print summary to get var explained
  summary(fitreg_total, standardized = T, rsquare = T)
  
  
  ## Total cholesterol EpiScore - regress total cholesterol EpiScore, age and sex on cog function intercept 

  
  reg_totalEpi <- '
 Ig ~ Total_chol_Epi + Age_w1 + sex
'
  
  fitreg_totalEpi <- growth(model = c(pgmodel1, pgmodel2, pgmodel3, pgmodel4, 
                                    pgmodel5, pgmodel6, pgmodel7, pgmodel8, 
                                    pgmodel9, pgmodel10, pgmodel11, pgmodel12, 
                                    pgmodel13, general_4p, reg_totalEpi), 
                          dset_mod,  missing = "ml.x")
  #print summary to get var explained
  summary(fitreg_totalEpi, standardized = T, rsquare = T)
  
  
  ## Total cholesterol and Total cholesterol EpiScore - 
  #regress total cholesterol, total cholesterol EpiScore, age and sex on cog function intercept 
  
  
  reg_total_totalEpi <- '
 Ig ~ bld_choles_w1 + Total_chol_Epi + Age_w1 + sex
'
  
  fitreg_total_totalEpi <- growth(model = c(pgmodel1, pgmodel2, pgmodel3, 
                                        pgmodel4, pgmodel5, pgmodel6, 
                                        pgmodel7, pgmodel8, pgmodel9, pgmodel10, 
                                        pgmodel11, pgmodel12, pgmodel13, general_4p, 
                                        reg_total_totalEpi), 
                              dset_mod,  missing = "ml.x")
  #print summary to get var explained
  summary(fitreg_total_totalEpi, standardized = T, rsquare = T)
  
  
  
  
  
  
  
