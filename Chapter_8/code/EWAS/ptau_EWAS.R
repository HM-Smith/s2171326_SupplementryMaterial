#required libraries
library(dplyr)
library(lme4)
library(lmerTest)
library(tictoc)
library(foreach)
library(doParallel)

#Load methylation data
load("20180128_0046_betas_ct_batch.Rdata")

#load phenotype file 
load("20191125_1614_pheno_combined.Rdata")


#sample names should be rows and CpGs columns so data needs transposing
meth_data = t(betas_ct_batch)

#remove harmony and parkinsons studies
pheno_combined <- pheno_combined[!(pheno_combined$STUDY %in% c("Harmony", "Parkinson")), ]

#read in ptau217, rename columns and select the ones of interest
response = read.csv("pTau217.csv")
names(response)[2] <- "TWINNR"
names(response)[6] <- "STUDY"

response = response[, c(2,6,9)]

response$TWINNR = as.character(response$TWINNR)

#remove studies without both methylation and proteins
response = response[!(response$STUDY %in% c("SATSA7")), ]
pheno_combined <- pheno_combined[!(pheno_combined$STUDY %in% c("IPT3")), ]


#rename studies
response$STUDY <- recode(response$STUDY,
                         "SATSA10" = "IPT10",
                         "SATSA5"  = "IPT5",
                         "SATSA6"  = "IPT6",
                         "SATSA8"  = "IPT8",
                         "SATSA9"  = "IPT9")

#combine datasets
response = left_join(pheno_combined, response, by = c("TWINNR", "STUDY"))


#subset to complete data
response = response[complete.cases(response), ]

#set rownames as sample names
rownames(response) = response$Sample

#merge with methylation data
meth_data = merge(response, meth_data, by=0, all.x=TRUE)


#create list of CpGs
CpGs <- names(meth_data)[14:ncol(meth_data)]

cores <- detectCores()
cl <- makeCluster(30)
registerDoParallel(cl)


# Set up the results list to collect the results
results_list <- foreach(j = CpGs, .packages = c("lmerTest", "lme4"), .combine = 'rbind') %dopar% {
  # Linear regression model
  mod <- lmer(scale(log_pTau217) ~ scale(meth_data[,j]) + 
                scale(AGE) + 
                as.factor(SEX) +
                (1|PAIRID) +
                (1|TWINNR) +
                as.factor(CHIP) + 
                as.factor(SMOKE), data = meth_data)
  
  mod.sum <- summary(mod)$coefficients
  mod.con <- confint(mod)
  
  # Create a result for each CpG
  data.frame(
    CpG = j,
    Protein = "pTau217",  # Modify as needed
    N = nobs(mod),
    Beta = mod.sum["scale(meth_data[, j])", 1],
    SE = mod.sum["scale(meth_data[, j])", 2],
    P.value = mod.sum["scale(meth_data[, j])", 5],
    Lower.CI = mod.con["scale(meth_data[, j])", 1],
    Upper.CI = mod.con["scale(meth_data[, j])", 2]
  )
}

#save out data
write.csv(results_list, file = "pTau217_EWAS_TWINNR.csv")

q()