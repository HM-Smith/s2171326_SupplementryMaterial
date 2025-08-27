# Load required libraries
library(data.table)
library(dplyr)

# Load target file with cohort/sample information
d = readRDS("/targets_3489_bloodonly.rds")

# Subset to LBC1936 cohort, wave 2, set 1
d36 <- d[d$cohort=="LBC36" & d$WAVE==2,]

# Load methylation data (betas) and transpose so rows = samples, cols = CpGs
dat <- readRDS("LBC_betas_3489_bloodonly.rds")
meth = t(dat)

# Subset methylation data to LBC1936 wave 2 samples
tmp36 <- which(rownames(meth) %in% d36$Basename)
meth36 <- meth[tmp36,] 

# Load in pre-computed CpG weights (not provided here, file path missing)
g_wts <- readRDS("") 

# Reformat weights: keep CpG IDs and corresponding betas
meanBetas <- g_wts[,c(2,1)]
names(meanBetas)=c("CpG","Beta")

# Convert methylation data to dataframe
meth36 = as.data.frame(meth36)

# Subset methylation data to CpGs present in weights file
a36 = which(names(meth36) %in% meanBetas$CpG) 
meth36 <- meth36[,a36] 

# Check for missing values (should return only FALSE if no NA’s)
table(is.na(meth36))
# FALSE
# 7920

# Convert methylation data back to a matrix
meth36 <- as.matrix(meth36) 

# Ensure CpGs in methylation data align with weights
b36 = which(meanBetas$CpG %in% colnames(meth36))
mean_betas <- meanBetas[b36,] 

# Order methylation data columns to match CpG order in weights
meth36 <- meth36[,match(colnames(meth36), mean_betas$CpG)]
meth36 <- meth36[,mean_betas$CpG]

# Convert back to dataframe for easy handling
meth36 = as.data.frame(meth36)

# Confirm CpG names match perfectly
identical(names(meth36), mean_betas$CpG) 
# [1] TRUE

# Convert again to matrix for matrix multiplication
meth36 = as.matrix(meth36)

# Create predictor score (linear combination of betas and weights)
pred36 <- meth36 %*% mean_betas$Beta
pred36 <- as.data.frame(pred36)
names(pred36) <- "GFAP_EpiScore"         # Rename column
pred36$Basename <- rownames(pred36)      # Add sample IDs

# Merge predictor score back into phenotype data
d36 <- full_join(d36, pred36, by="Basename")
d36 <- d36 %>% rename("lbc36no" = ID)

# Convert technical covariates to factors
d36$set <- as.factor(d36$set)
d36$date <- as.factor(d36$date)
d36$array <- as.factor(d36$array)

# Drop unused factor levels for "date"
d36$date <- droplevels(d36$date)

# Load lme4 for mixed effects modeling
library(lme4)

# Model 1: adjust EpiScore for white blood cell counts + technical effects
mod = lmer(GFAP_EpiScore ~ scale(neut) + scale(lymph) + scale(mono)+ scale(eosin)+
             scale(baso) +(1|set)  +(1|array) + (1|date), 
           na.action = na.exclude, data = d36)

# Save residuals = WBC + technical-adjusted score
d36$EpiScore_WBC_tech_adj = resid(mod)

# Model 2: adjust only for technical effects (no WBC adjustment)
mod2 = lmer(EpiScore ~ (1|set)  +(1|array) + (1|date), 
            na.action = na.exclude, data = d36)

# Save residuals = technical-adjusted score
d36$GFAP_EpiScore_tech_adj = resid(mod2)

# Write output dataset to CSV
write.csv(d36, file = "EpiScore.csv")
