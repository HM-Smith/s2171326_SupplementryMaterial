#### library ####
library("tidyverse")
library("lavaan")
library("foreach")
library("doParallel")
library("readr")

## read in LBC1936 data
data <- read.csv("LBC1921_WAVE1_cohort_data_KORA_EpiScore_Adjusted_31012023.csv")

#set sex as factor 
data$sex <- as.factor(data$sex)

#rename data from to work in script
dset_mod <- data

#scale age 
dset_mod$age <- as.vector(scale(dset_mod$age))


#### growth curve models for individual cognitive tests ####

pgmodel1 <- '
Ivftot =~ 1*vftot79 + 1*vftot83 + 1*vtotal87 + 1*vftotal_w4 + 1*vftot_w5
Svftot =~ 0*vftot79 + 4.3*vftot83 + 7.5*vtotal87 + 11*vftotal_w4 + 13*vftot_w5
'
fit1 <- growth(pgmodel1, dset_mod, missing = "ml.x")
#summary(fit1, standardized = T)

pgmodel2 <- '
Inart =~ 1*nart79 + 1*nart83 + 1*nart87 + 1*nart_w4
Snart =~ 0*nart79 + 4.3*nart83 + 7.5*nart87 + 11*nart_w4 
'
fit2 <- growth(pgmodel2, dset_mod, missing = "ml.x") #Snart has negative latent variance 
#summary(fit2, standardized = T)

pgmodel3 <- '
Ilogmem =~ 1*logmem79 + 1*logtot83 + 1*logtot87 + 1*logmem_w4 + 1*logmem_w5
Slogmem =~ 0*logmem79 + 4.3*logtot83 + 7.5*logtot87 + 11*logmem_w4 + 13*logmem_w5
'
fit3 <- growth(pgmodel3, dset_mod, missing = "ml.x")
#summary(fit3, standardized = T)

pgmodel4 <- '
Iravens =~ 1*ravens79 + 1*ravens83 + 1*ravens87 + 1*ravens_w4 + 1*ravens_w5
Sravens =~ 0*ravens79 + 4.3*ravens83 + 7.5*ravens87 + 11*ravens_w4 + 13*ravens_w5
'
fit4 <- growth(pgmodel4, dset_mod, missing = "ml.x")
#summary(fit4, standardized = T)


#create measurement model for general cognitive function and change
general_4p <- '

#latent variables 
Ig =~ 1*Ivftot + 0.889*Ilogmem + 0.914*Inart + 0.794*Iravens

Sg =~ 1*Svftot + 2.078*Slogmem + 0.177*Snart + 0.951*Sravens

#indicator as scaling reference: loading=1, int=0
Ivftot ~ 0*1
Svftot ~ 0*1 

#covariances 
Snart~~0*Snart
Sravens ~~ 0*Sravens
'

#fit model
fitGen_4p <- growth(model = c(pgmodel1,pgmodel2, pgmodel3,pgmodel4, general_4p), dset_mod,  missing = "ml.x", em.h1.iter.max = 1000)

#extract fitmeasures
fitmeasures(fitGen_4p, c("cfi", "tli", "RMSEA", "SRMR"))
#cfi   tli rmsea  srmr
#0.969 0.970 0.037 0.082


#model summary 
summary(fitGen_4p, standardized = T)


################################################################
#register 84 cores to allow parallel processing 
cores <- detectCores()
cl <- makeCluster(84)
registerDoParallel(cl)

#make list of EpiScores
list_e <- colnames(dset_mod[105:188])
episcores = colnames(dset_mod)[which(colnames(dset_mod)%in% list_e)]

#for loop that loops over episcores using parallel processing to perform regression in SEM
#saves out results for each episcore: G intercept ~ episcore (tmp) + age + sex
foreach(i = episcores, .packages = c("lavaan")) %dopar% { 
  
  dset_mod$tmp = dset_mod[,i]
  
  reg_Ig <- '
 Ig ~ tmp + age + sex
'
  
  fitreg_g <- growth(model = c(pgmodel1, pgmodel2, pgmodel3, pgmodel4, general_4p, reg_Ig), dset_mod,  missing = "ml.x")
  
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
  
  write.csv(results, paste0(i, "_G_intercept_No_Domain_LBC1921_assocs_84_EpiScores_basic_21042023.csv"), row.names = F)
  
  
}

##########################################################################################################################################################################################################################################################
##########################################################################################################################################################################################################################################################

#make list of EpiScore
list_e <- colnames(dset_mod[105:188])
episcores = colnames(dset_mod)[which(colnames(dset_mod)%in% list_e)]


#for loop that loops over episcores using parallel processing to perform regression in SEM
#saves out results for each episcore: G slope ~ episcore (tmp) + age + sex
foreach(i = episcores, .packages = c("lavaan")) %dopar% { 
  
  dset_mod$tmp = dset_mod[,i]
  
  reg_Sg <- '
 Sg ~ tmp + age + sex
'
  
  fitreg_g <- growth(model = c(pgmodel1, pgmodel2, pgmodel3, pgmodel4, general_4p, reg_Sg), dset_mod,  missing = "ml.x")
  
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
  
  write.csv(results, paste0(i, "_G_slope_No_Domain_LBC1921_assocs_84_EpiScores_basic_21042023.csv"), row.names = F)
  
  print(i)
}


