#### library ####
library("tidyverse")
library("lavaan")
library("foreach")
library("doParallel")
library("readr")

#read in congitive/episcore data 
data <- read.csv("GS20K_Cognitive_KORA_adjust_EpiScores_combined_20012023.csv")

#set sex as factor
data$sex <- as.factor(data$sex)

#create measurement model of general cognitive function
G_model <- '

 #latent variables
 G =~ LM + verbal_total + digit_symbol + vocabulary


 '

#fit model 
fitGen_4p <- sem(model = G_model, data,  missing = "ml.x")

#extract fitmeausres
fitmeasures(fitGen_4p, c("cfi", "tli", "RMSEA", "SRMR")) #cfi   tli rmsea  srmr 0.734 0.202 0.219 0.066

#model summary
summary(fitGen_4p, standardized = T)

################################################################
#register 84 cores to allow parallel processing 
cores <- detectCores()
cl <- makeCluster(84)
registerDoParallel(cl)


#make list of episcores
list_e <- colnames(data[25:108])
episcores = colnames(data)[which(colnames(data)%in% list_e)]


#for loop that loops over episcores using parallel processing to perform regression in SEM
#saves out results for each episcore: G ~ episcore (tmp) + age + sex
foreach(i = episcores, .packages = c("lavaan")) %dopar% { 
  
  data$tmp = data[,i]
  
  reg_Ig <- '
 G ~ tmp + age + sex
'
  
  fitreg_g <- sem(model = c(G_model, reg_Ig), data,  missing = "ml.x")
  
  output = standardizedSolution(fitreg_g)
  ind = which(output$lhs == "G" & output$op == "~" & output$rhs=="tmp")
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
  
  write.csv(results, paste0(i, "_G__Generation_Scot_assocs_84_EpiScores_basic_20012023.csv"))
  
  
}
