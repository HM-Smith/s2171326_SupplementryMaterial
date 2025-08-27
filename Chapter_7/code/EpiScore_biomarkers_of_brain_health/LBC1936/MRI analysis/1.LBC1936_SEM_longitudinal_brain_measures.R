library(readr)
library(lavaan)
library(tidyverse)
library(dplyr)

#read in MRI data
LBC <- read.csv("LBC1936_WAVE1_cohort_data_KORA_EpiScore_Adjusted_17072023.csv")

#set sex as factor 
LBC$sex <- as.factor(LBC$sex)


## calculating time between years 
data2 <- LBC %>% filter(Age_w5 != "NA")
age_w1.1 <- colMeans(data2["Age_w1"], na.rm = T)
age_w2.1 <- colMeans(data2["Age_w2"], na.rm = T)
age_w3.1 <- colMeans(data2["Age_w3"], na.rm = T)
age_w4.1 <- colMeans(data2["Age_w4"], na.rm = T)
age_w5.1 <- colMeans(data2["Age_w5"], na.rm = T)

lag1 <- age_w3.1 - age_w2.1
lag2 <- age_w4.1 - age_w2.1
lag3 <- age_w5.1 - age_w2.1


print(lag1)#3.76
print(lag2)#6.83
print(lag3)#9.55


#scale age and intracrainal volume 
LBC$ICV_mm3_wX = as.vector(scale(LBC$ICV_mm3_wX))
LBC$Age_w1 = as.vector(scale(LBC$Age_w1))
LBC$Age_w2 = as.vector(scale(LBC$Age_w2))
LBC$Age_w3 = as.vector(scale(LBC$Age_w3))
LBC$Age_w4 = as.vector(scale(LBC$Age_w4))
LBC$Age_w5 = as.vector(scale(LBC$Age_w5))


#rescale MRI variables to aid in model convergence 
LBC$gm_mm3_w2 <- LBC$gm_mm3_w2/10000
LBC$gm_mm3_w3 <- LBC$gm_mm3_w3/10000
LBC$gm_mm3_w4 <- LBC$gm_mm3_w4/10000
LBC$gm_mm3_w5 <- LBC$gm_mm3_w5/10000

LBC$brain_mm3_w2 <- LBC$brain_mm3_w2/10000
LBC$brain_mm3_w3 <- LBC$brain_mm3_w3/10000
LBC$brain_mm3_w4 <- LBC$brain_mm3_w4/10000
LBC$brain_mm3_w5 <- LBC$brain_mm3_w5/10000

LBC$nawm_mm3_w2 <- LBC$nawm_mm3_w2/10000
LBC$nawm_mm3_w3 <- LBC$nawm_mm3_w3/10000
LBC$nawm_mm3_w4 <- LBC$nawm_mm3_w4/10000
LBC$nawm_mm3_w5 <- LBC$nawm_mm3_w5/10000

attach(LBC)
LBC$wmh_mm3_log_w2 <- log(wmh_mm3_w2 + 1)
LBC$wmh_mm3_log_w3 <- log(wmh_mm3_w3 + 1)
LBC$wmh_mm3_log_w4 <- log(wmh_mm3_w4 + 1)
LBC$wmh_mm3_log_w5 <- log(wmh_mm3_w5 + 1)
detach()


########################################################################################
########################################################################################
#create grey matter measurement model for intercept and slope 
  GM_model <- '
   # latent variables
     IGMV =~ 1*gm_mm3_w2 + 1*gm_mm3_w3 + 1*gm_mm3_w4 + 1*gm_mm3_w5
     SGMV =~ 0*gm_mm3_w2 + 3.76*gm_mm3_w3 + 6.83*gm_mm3_w4 + 9.55*gm_mm3_w5
  
   # regressions 
     gm_mm3_w2 ~ Age_w2
     gm_mm3_w3 ~ Age_w3 
     gm_mm3_w4 ~ Age_w4 
     gm_mm3_w5 ~ Age_w5 
    
'

#fit model 
fit_GM <- growth(model = GM_model, data  = LBC, missing = "ml.x")

#model summary 
summary(fit_GM, standardized = TRUE, fit.measures = TRUE)

#extract fit measures 
fitmeasures(fit_GM, c("cfi", "tli", "RMSEA", "SRMR"))
#cfi   tli rmsea  srmr
#0.941 0.924 0.101 0.035


#make list of episcores 
list_e <- colnames(LBC)[335:418]
episcores = colnames(LBC)[which(colnames(LBC)%in% list_e)]

#create two empty results dataframes 
results <- data.frame(SeqId = list_e, n = NA, beta = NA , SE = NA, P = NA, ci.upper = NA, ci.lower = NA, cfi = NA, rmsea = NA, srmr = NA, tli = NA, phenotype = NA)
results_slope <- data.frame(SeqId = list_e, n = NA, beta = NA , SE = NA, P = NA, ci.upper = NA, ci.lower = NA, cfi = NA, rmsea = NA, srmr = NA, tli = NA, phenotype = NA)

#set episcores as rownames 
rownames(results) = results$SeqId # This will allow you to index with results[i,]
rownames(results_slope) = results_slope$SeqId

#for loop that loops over episcores and performs regression analysis
for(i in episcores) { 
  LBC$tmp = LBC[,i]


    GM_model_reg <-  '
    # latent variables
    IGMV =~ 1*gm_mm3_w2 + 1*gm_mm3_w3 + 1*gm_mm3_w4 + 1*gm_mm3_w5
    SGMV =~ 0*gm_mm3_w2 + 3.76*gm_mm3_w3 + 6.83*gm_mm3_w4 + 9.55*gm_mm3_w5
    
    # regressions 
    gm_mm3_w2 ~ Age_w2 
    gm_mm3_w3 ~ Age_w3 
    gm_mm3_w4 ~ Age_w4 
    gm_mm3_w5 ~ Age_w5 
    tmp ~ Age_w1
    IGMV ~ tmp + sex + ICV_mm3_wX
    SGMV ~ tmp + sex 
    
    '
             
    GM_fit_reg <- growth(model = GM_model_reg, data = LBC, missing = "ml.x" )


  output = standardizedSolution(GM_fit_reg)
  ind = which(output$lhs == "IGMV" & output$op == "~" & output$rhs=="tmp")
  snd = which(output$lhs == "SGMV" & output$op == "~" & output$rhs=="tmp")
  
  Beta = output[ind, "est.std"] 
  SE <- output[ind, "se"]
  p <- output[ind, "pvalue"]
  n <- nobs(GM_fit_reg)
  ci.upper <- output[ind, "ci.upper"]
  ci.lower <- output[ind, "ci.lower"]
  fitmeasures <- fitMeasures(GM_fit_reg, c("cfi","rmsea","srmr", "tli"))
  cfi <- fitmeasures["cfi"]
  rmsea <- fitmeasures["rmsea"]
  srmr <- fitmeasures["srmr"]
  tli <- fitmeasures["tli"]
  pheno <- output[ind, "lhs"]
  
  results[i,1] <- i
  results[i,2] <- n
  results[i,3] <- Beta
  results[i,4] <- SE
  results[i,5] <- p
  results[i,6] <- ci.upper
  results[i,7] <- ci.lower
  results[i,8] <- cfi
  results[i,9] <- rmsea
  results[i,10] <- srmr
  results[i,11] <- tli
  results[i,12] <- pheno
  
  Beta_slope= output[snd, "est.std"] 
  SE_slope <- output[snd, "se"]
  p_slope <- output[snd, "pvalue"]
  n_slope <- nobs(GM_fit_reg)
  ci.upper_slope <- output[snd, "ci.upper"]
  ci.lower_slope <- output[snd, "ci.lower"]
  fitmeasures_slope <- fitMeasures(GM_fit_reg, c("cfi","rmsea","srmr", "tli"))
  cfi_slope <- fitmeasures_slope["cfi"]
  rmsea_slope <- fitmeasures_slope["rmsea"]
  srmr_slope <- fitmeasures_slope["srmr"]
  tli_slope <- fitmeasures_slope["tli"]
  pheno_slope <- output[snd, "lhs"]
  
  results_slope[i,1] <- i
  results_slope[i,2] <- n_slope
  results_slope[i,3] <- Beta_slope
  results_slope[i,4] <- SE_slope
  results_slope[i,5] <- p_slope
  results_slope[i,6] <- ci.upper_slope
  results_slope[i,7] <- ci.lower_slope
  results_slope[i,8] <- cfi_slope
  results_slope[i,9] <- rmsea_slope
  results_slope[i,10] <- srmr_slope
  results_slope[i,11] <- tli_slope
  results_slope[i, 12] <- pheno_slope
  
  print(i)
  
  
  }



#save out results
write.csv(results, file = "GreyMatter_intercept_LBC1936_assocs_84_EpiScores_basic_15092023.csv")
write.csv(results_slope, file = "GreyMatter_slope_LBC1936_assocs_84_EpiScores_basic_15092023.csv")

#remove results dataframes
rm(results)
rm(results_slope)


###########################################################################################
###########################################################################################
#create total brain volume measurement model for intercept and slope 
Brain_model <- '
   # latent variables
     ITMV =~ 1*brain_mm3_w2 + 1*brain_mm3_w3 + 1*brain_mm3_w4 + 1*brain_mm3_w5
     STMV =~ 0*brain_mm3_w2 + 3.76*brain_mm3_w3 + 6.83*brain_mm3_w4 + 9.55*brain_mm3_w5
  
   # regressions 
     brain_mm3_w2 ~ Age_w2 
     brain_mm3_w3 ~ Age_w3 
     brain_mm3_w4 ~ Age_w4 
     brain_mm3_w5 ~ Age_w5 
    
'

#fit model
Brain_fit <- growth(model = Brain_model, data  = LBC, missing = "ml.x")

#model summary 
summary(Brain_fit, standardized = TRUE, fit.measures = TRUE)


#extract fit measures 
fitmeasures(Brain_fit, c("cfi", "tli", "RMSEA", "SRMR"))
#cfi   tli rmsea  srmr
#0.992 0.990 0.044 0.016


#make list of episcores
list_e <- colnames(LBC)[335:418]
episcores = colnames(LBC)[which(colnames(LBC)%in% list_e)]

#make two empty results dataframe 
results <- data.frame(SeqId = list_e, n = NA, beta = NA , SE = NA, P = NA, ci.upper = NA, ci.lower = NA, cfi = NA, rmsea = NA, srmr = NA, tli = NA, phenotype = NA)
results_slope <- data.frame(SeqId = list_e, n = NA, beta = NA , SE = NA, P = NA, ci.upper = NA, ci.lower = NA, cfi = NA, rmsea = NA, srmr = NA, tli = NA, phenotype = NA)

#set episcores as rownames 
rownames(results) = results$SeqId # This will allow you to index with results[i,]
rownames(results_slope) = results_slope$SeqId


#for loop that loops over episcores and performs regression analysis
for(i in episcores) { 
  LBC$tmp = LBC[,i]
  
  
  Brain_model_reg <-  '
    # latent variables
     ITMV =~ 1*brain_mm3_w2 + 1*brain_mm3_w3 + 1*brain_mm3_w4 + 1*brain_mm3_w5
     STMV =~ 0*brain_mm3_w2 + 3.76*brain_mm3_w3 + 6.83*brain_mm3_w4 + 9.55*brain_mm3_w5
  
   # regressions 
     brain_mm3_w2 ~ Age_w2 
     brain_mm3_w3 ~ Age_w3 
     brain_mm3_w4 ~ Age_w4 
     brain_mm3_w5 ~ Age_w5 
     tmp ~ Age_w1
     ITMV ~ tmp + sex + ICV_mm3_wX
     STMV ~ tmp + sex
    '
  
  Brain_fit_reg <- growth(model = Brain_model_reg, data = LBC, missing = "ml.x" )
  
  
  output = standardizedSolution(Brain_fit_reg)
  ind = which(output$lhs == "ITMV" & output$op == "~" & output$rhs=="tmp")
  snd = which(output$lhs == "STMV" & output$op == "~" & output$rhs=="tmp")
  
  Beta = output[ind, "est.std"] 
  SE <- output[ind, "se"]
  p <- output[ind, "pvalue"]
  n <- nobs(Brain_fit_reg)
  ci.upper <- output[ind, "ci.upper"]
  ci.lower <- output[ind, "ci.lower"]
  fitmeasures <- fitMeasures(Brain_fit_reg, c("cfi","rmsea","srmr", "tli"))
  cfi <- fitmeasures["cfi"]
  rmsea <- fitmeasures["rmsea"]
  srmr <- fitmeasures["srmr"]
  tli <- fitmeasures["tli"]
  pheno <- output[ind, "lhs"]
  
  results[i,1] <- i
  results[i,2] <- n
  results[i,3] <- Beta
  results[i,4] <- SE
  results[i,5] <- p
  results[i,6] <- ci.upper
  results[i,7] <- ci.lower
  results[i,8] <- cfi
  results[i,9] <- rmsea
  results[i,10] <- srmr
  results[i,11] <- tli
  results[i,12] <- pheno
  
  Beta_slope= output[snd, "est.std"] 
  SE_slope <- output[snd, "se"]
  p_slope <- output[snd, "pvalue"]
  n_slope <- nobs(Brain_fit_reg)
  ci.upper_slope <- output[snd, "ci.upper"]
  ci.lower_slope <- output[snd, "ci.lower"]
  fitmeasures_slope <- fitMeasures(Brain_fit_reg, c("cfi","rmsea","srmr", "tli"))
  cfi_slope <- fitmeasures_slope["cfi"]
  rmsea_slope <- fitmeasures_slope["rmsea"]
  srmr_slope <- fitmeasures_slope["srmr"]
  tli_slope <- fitmeasures_slope["tli"]#
  pheno_slope <- output[snd, "lhs"]
  
  results_slope[i,1] <- i
  results_slope[i,2] <- n_slope
  results_slope[i,3] <- Beta_slope
  results_slope[i,4] <- SE_slope
  results_slope[i,5] <- p_slope
  results_slope[i,6] <- ci.upper_slope
  results_slope[i,7] <- ci.lower_slope
  results_slope[i,8] <- cfi_slope
  results_slope[i,9] <- rmsea_slope
  results_slope[i,10] <- srmr_slope
  results_slope[i,11] <- tli_slope
  results_slope[i,12] <- pheno_slope
  
  print(i)
  
}

#save out results
write.csv(results, file = "Totalbrain_intercept_LBC1936_assocs_84_EpiScores_basic_15092023.csv")
write.csv(results_slope, file ="Totalbrain_slope_LBC1936_assocs_84_EpiScores_basic_15092023.csv")

#remove results dataframe
rm(results)
rm(results_slope)

###########################################################################################
###########################################################################################
#create normal appearing white matter volume measurement model for intercept and slope 
NAWM_model <- '
   # latent variables
     INAWMV =~ 1*nawm_mm3_w2 + 1*nawm_mm3_w3 + 1*nawm_mm3_w4 + 1*nawm_mm3_w5
     SNAWMV =~ 0*nawm_mm3_w2 + 3.76*nawm_mm3_w3 + 6.83*nawm_mm3_w4 + 9.55*nawm_mm3_w5
  
   # regressions 
     nawm_mm3_w2 ~ Age_w2 
     nawm_mm3_w3 ~ Age_w3 
     nawm_mm3_w4 ~ Age_w4 
     nawm_mm3_w5 ~ Age_w5 
    
'

#fit model
fit_NAWM <- growth(model = NAWM_model, data  = LBC, missing = "ml.x")

#model summary
summary(fit_NAWM, standardized = TRUE, fit.measures = TRUE)

#extract fit measures
fitmeasures(fit_NAWM, c("cfi", "tli", "RMSEA", "SRMR"))
#cfi   tli rmsea  srmr
#0.994 0.992 0.034 0.017


#make list of episcores 
list_e <- colnames(LBC)[335:418]
episcores = colnames(LBC)[which(colnames(LBC)%in% list_e)]

#make two empty results dataframe
results <- data.frame(SeqId = list_e, n = NA, beta = NA , SE = NA, P = NA, ci.upper = NA, ci.lower = NA, cfi = NA, rmsea = NA, srmr = NA, tli = NA, phenotype = NA)
results_slope <- data.frame(SeqId = list_e, n = NA, beta = NA , SE = NA, P = NA, ci.upper = NA, ci.lower = NA, cfi = NA, rmsea = NA, srmr = NA, tli = NA, phenotype = NA)

#set episcores as rownames 
rownames(results) = results$SeqId # This will allow you to index with results[i,]
rownames(results_slope) = results_slope$SeqId


#for loop that loops over episcores and performs regression analysis
for(i in episcores) { 
  LBC$tmp = LBC[,i]
  
  
  NAWM_model_reg <-  '
  # latent variables
  INAWMV =~ 1*nawm_mm3_w2 + 1*nawm_mm3_w3 + 1*nawm_mm3_w4 + 1*nawm_mm3_w5
  SNAWMV =~ 0*nawm_mm3_w2 + 3.76*nawm_mm3_w3 + 6.83*nawm_mm3_w4 + 9.55*nawm_mm3_w5
  
  # regressions 
  nawm_mm3_w2 ~ Age_w2 
  nawm_mm3_w3 ~ Age_w3 
  nawm_mm3_w4 ~ Age_w4 
  nawm_mm3_w5 ~ Age_w5 
  tmp ~ Age_w1
  INAWMV ~ tmp + sex + ICV_mm3_wX
  SNAWMV ~ tmp + sex
  '
  
  NAWM_fit_reg <- growth(model = NAWM_model_reg, data = LBC, missing = "ml.x" )
  
  
  output = standardizedSolution(NAWM_fit_reg)
  ind = which(output$lhs == "INAWMV" & output$op == "~" & output$rhs=="tmp")
  snd = which(output$lhs == "SNAWMV" & output$op == "~" & output$rhs=="tmp")
  
  Beta = output[ind, "est.std"] 
  SE <- output[ind, "se"]
  p <- output[ind, "pvalue"]
  n <- nobs(NAWM_fit_reg)
  ci.upper <- output[ind, "ci.upper"]
  ci.lower <- output[ind, "ci.lower"]
  fitmeasures <- fitMeasures(NAWM_fit_reg, c("cfi","rmsea","srmr", "tli"))
  cfi <- fitmeasures["cfi"]
  rmsea <- fitmeasures["rmsea"]
  srmr <- fitmeasures["srmr"]
  tli <- fitmeasures["tli"]
  pheno <- output[ind, "lhs"]
  
  results[i,1] <- i
  results[i,2] <- n
  results[i,3] <- Beta
  results[i,4] <- SE
  results[i,5] <- p
  results[i,6] <- ci.upper
  results[i,7] <- ci.lower
  results[i,8] <- cfi
  results[i,9] <- rmsea
  results[i,10] <- srmr
  results[i,11] <- tli
  results[i,12] <- pheno
  
  Beta_slope= output[snd, "est.std"] 
  SE_slope <- output[snd, "se"]
  p_slope <- output[snd, "pvalue"]
  n_slope <- nobs(NAWM_fit_reg)
  ci.upper_slope <- output[snd, "ci.upper"]
  ci.lower_slope <- output[snd, "ci.lower"]
  fitmeasures_slope <- fitMeasures(NAWM_fit_reg, c("cfi","rmsea","srmr", "tli"))
  cfi_slope <- fitmeasures_slope["cfi"]
  rmsea_slope <- fitmeasures_slope["rmsea"]
  srmr_slope <- fitmeasures_slope["srmr"]
  tli_slope <- fitmeasures_slope["tli"]
  pheno_slope <- output[snd, "lhs"]
  
  results_slope[i,1] <- i
  results_slope[i,2] <- n_slope
  results_slope[i,3] <- Beta_slope
  results_slope[i,4] <- SE_slope
  results_slope[i,5] <- p_slope
  results_slope[i,6] <- ci.upper_slope
  results_slope[i,7] <- ci.lower_slope
  results_slope[i,8] <- cfi_slope
  results_slope[i,9] <- rmsea_slope
  results_slope[i,10] <- srmr_slope
  results_slope[i,11] <- tli_slope
  results_slope[i,12] <- pheno_slope
  
  print(i)
}

#save out results 
write.csv(results, file = "NAWM_intercept_LBC1936_assocs_84_EpiScores_basic_15092023.csv")
write.csv(results_slope, file ="NAWM_slope_LBC1936_assocs_84_EpiScores_basic_15092023.csv")

#remove results dataframes
rm(results)
rm(results_slope)

###########################################################################################
###########################################################################################
#create white matter hyperintensity volume measurement model for intercept and slope 
WMH_model <- '
   # latent variables
     IWMHV =~ 1*wmh_mm3_log_w2 + 1*wmh_mm3_log_w3 + 1*wmh_mm3_log_w4 + 1*wmh_mm3_log_w5
     SWMHV =~ 0*wmh_mm3_log_w2 + 3.76*wmh_mm3_log_w3 + 6.83*wmh_mm3_log_w4 + 9.55*wmh_mm3_log_w5
  
   # regressions 
     wmh_mm3_log_w2 ~ Age_w2 
     wmh_mm3_log_w3 ~ Age_w3 
     wmh_mm3_log_w4 ~ Age_w4
     wmh_mm3_log_w5 ~ Age_w5 
    
'

#fit model
WMH_fit <- growth(model = WMH_model, data  = LBC, missing = "ml.x")

#model summary
summary(WMH_fit, standardized = TRUE, fit.measures = TRUE)

#extract fit measures
fitmeasures(WMH_fit, c("cfi", "tli", "RMSEA", "SRMR"))
#cfi   tli rmsea  srmr
#0.937 0.919 0.099 0.062


#make list of episcores
list_e <- colnames(LBC)[335:418]
episcores = colnames(LBC)[which(colnames(LBC)%in% list_e)]

#make two empty results dataframes
results <- data.frame(SeqId = list_e, n = NA, beta = NA , SE = NA, P = NA, ci.upper = NA, ci.lower = NA, cfi = NA, rmsea = NA, srmr = NA, tli = NA, phenotype = NA)
results_slope <- data.frame(SeqId = list_e, n = NA, beta = NA , SE = NA, P = NA, ci.upper = NA, ci.lower = NA, cfi = NA, rmsea = NA, srmr = NA, tli = NA, phenotype = NA)

#set episcores as rownames 
rownames(results) = results$SeqId # This will allow you to index with results[i,]
rownames(results_slope) = results_slope$SeqId


#for loop that loops over episcores and performs regression analysis
for(i in episcores) { 
  LBC$tmp = LBC[,i]
  
  
  WMH_model_reg <-  '
  # latent variables
     IWMHV =~ 1*wmh_mm3_log_w2 + 1*wmh_mm3_log_w3 + 1*wmh_mm3_log_w4 + 1*wmh_mm3_log_w5
     SWMHV =~ 0*wmh_mm3_log_w2 + 3.76*wmh_mm3_log_w3 + 6.83*wmh_mm3_log_w4 + 9.55*wmh_mm3_log_w5
  
   # regressions 
     wmh_mm3_log_w2 ~ Age_w2 
     wmh_mm3_log_w3 ~ Age_w3 
     wmh_mm3_log_w4 ~ Age_w4 
     wmh_mm3_log_w5 ~ Age_w5 
     tmp ~ Age_w1
     IWMHV ~ tmp + sex + ICV_mm3_wX
     SWMHV ~ tmp + sex
  '
  
  WMH_fit_reg <- growth(model = WMH_model_reg, data = LBC, missing = "ml.x" )
  
  
  output = standardizedSolution(WMH_fit_reg)
  ind = which(output$lhs == "IWMHV" & output$op == "~" & output$rhs=="tmp")
  snd = which(output$lhs == "SWMHV" & output$op == "~" & output$rhs=="tmp")
  
  Beta = output[ind, "est.std"] 
  SE <- output[ind, "se"]
  p <- output[ind, "pvalue"]
  n <- nobs(WMH_fit_reg)
  ci.upper <- output[ind, "ci.upper"]
  ci.lower <- output[ind, "ci.lower"]
  fitmeasures <- fitMeasures(WMH_fit_reg, c("cfi","rmsea","srmr", "tli"))
  cfi <- fitmeasures["cfi"]
  rmsea <- fitmeasures["rmsea"]
  srmr <- fitmeasures["srmr"]
  tli <- fitmeasures["tli"]
  pheno <- output[ind, "lhs"]
  
  results[i,1] <- i
  results[i,2] <- n
  results[i,3] <- Beta
  results[i,4] <- SE
  results[i,5] <- p
  results[i,6] <- ci.upper
  results[i,7] <- ci.lower
  results[i,8] <- cfi
  results[i,9] <- rmsea
  results[i,10] <- srmr
  results[i,11] <- tli
  results[i,12] <- pheno
  
  Beta_slope= output[snd, "est.std"] 
  SE_slope <- output[snd, "se"]
  p_slope <- output[snd, "pvalue"]
  n_slope <- nobs(WMH_fit_reg)
  ci.upper_slope <- output[snd, "ci.upper"]
  ci.lower_slope <- output[snd, "ci.lower"]
  fitmeasures_slope <- fitMeasures(WMH_fit_reg, c("cfi","rmsea","srmr", "tli"))
  cfi_slope <- fitmeasures_slope["cfi"]
  rmsea_slope <- fitmeasures_slope["rmsea"]
  srmr_slope <- fitmeasures_slope["srmr"]
  tli_slope <- fitmeasures_slope["tli"]
  pheno_slope <- output[snd, "lhs"]
  
  results_slope[i,1] <- i
  results_slope[i,2] <- n_slope
  results_slope[i,3] <- Beta_slope
  results_slope[i,4] <- SE_slope
  results_slope[i,5] <- p_slope
  results_slope[i,6] <- ci.upper_slope
  results_slope[i,7] <- ci.lower_slope
  results_slope[i,8] <- cfi_slope
  results_slope[i,9] <- rmsea_slope
  results_slope[i,10] <- srmr_slope
  results_slope[i,11] <- tli_slope
  results_slope[i, 12] <- pheno_slope
  
}

#save out results
write.csv(results, file = "WMH_intercept_LBC1936_assocs_84_EpiScores_basic_15092023.csv")
write.csv(results_slope, file ="WMH_slope_LBC1936_assocs_84_EpiScores_basic_15092023.csv")


