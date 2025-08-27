#### library ####
library("tidyverse")
library("lavaan")
library("foreach")
library("doParallel")
library("readr")

#read in cognitive/episcore data 
data <- read.csv("GS20K_Cognitive_KORA_adjust_EpiScores_combined_20012023.csv")

#set sex as factor
data$sex <- as.factor(data$sex)

#read in covariates files
alc <- read.csv("2023-02-01_alc.csv")
smid <- read.csv("2023-02-13_simd.csv")
BMI <- read.csv("body.csv")

#rename id column
names(alc)[1] <- "Sample_Name"
names(smid)[1] <- "Sample_Name"
names(BMI)[1] <- "Sample_Name"

#read in EpiSmokEr scores and combine
w1 <- readRDS("wave1_epismoker.rds")
w3 <- readRDS("wave3_epismoker.rds")
w4 <- readRDS("wave4_epismoker.rds")
bind <- rbind(w1, w3)
bind <- rbind(bind, w4) # 18779 individuals with DNAm epismoker calculated 
bind$Sample_Sentrix_ID <- row.names(bind)
bind <- bind[-1]

#join cognitive data with covariates
data <- data %>% left_join(alc, by ="Sample_Name")
data <- data %>% left_join(smid, by ="Sample_Name")
data <- data %>% left_join(BMI, by ="Sample_Name")
data <- data %>% left_join(bind, by ="Sample_Sentrix_ID")

#scale covariates
data$units <- scale(data$units)
data$age <- scale(data$units)
data$rank <- scale(data$rank)
data$bmi <- scale(data$bmi)
data$smokingScore <- scale(data$smokingScore)


#create measurement model of general cognitive function
G_model <- '

 #latent variables
 G =~ LM + verbal_total + digit_symbol + vocabulary


 '

#run model 
fitGen_4p <- sem(model = G_model, data,  missing = "ml.x") 

#extract fitmeasures
fitmeasures(fitGen_4p, c("cfi", "tli", "RMSEA", "SRMR")) #cfi   tli rmsea  srmr 0.734 0.202 0.219 0.066

#summarise model 
summary(fitGen_4p, standardized = T)
################################################################
#register 84 cores for parallel processing 
cores <- detectCores()
cl <- makeCluster(84)
registerDoParallel(cl)

#make list of episcores
list_e <- colnames(data[25:108])
episcores = colnames(data)[which(colnames(data)%in% list_e)]

#for loop that loops over episcores using parallel processing to perform regression in SEM
#saves out results for each episcore: G ~ episcore (tmp) + age + sex + BMI + deprivation index + alcohol units + epismoker 
foreach(i = episcores, .packages = c("lavaan")) %dopar% { 
  
  data$tmp = data[,i]
  
  reg_Ig <- '
 G ~ tmp + age + sex + units + rank + bmi + smokingScore 
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
  
  write.csv(results, paste0(i, "_G__Generation_Scot_assocs_84_EpiScores_fullmodel_no_yrsedu_13022023.csv"))
  
  
}

