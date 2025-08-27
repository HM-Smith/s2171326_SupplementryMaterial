# Set working directory 
setwd("/Cluster_Filespace/Marioni_Group/GS/GS_methylation/GS20k/Chromosomes/") 

# Load requisite libraries
library(data.table)
library(limma)
library(coxme)
library(kinship2)

# Create function for mean imputation (not ideal but will do as sensitivity analysis for practical/computation reasons)
meanimpute <- function(x) ifelse(is.na(x),mean(x,na.rm=T),x)

# Create function for Outlier removal - 4SD
outlierID <- function(x, cut=4) {
  xx <- scale(x)
  retval <- ifelse (abs(xx) >= cut, TRUE, FALSE)
  retval
}

outlierTrim <- function(x, cut=4) {
  id <- outlierID(x, cut=cut)
  retval <- ifelse(id, NA, x)
}


#########################################################
######## STEP 1 - PREPARATION OF PHENOTYPE FILES ########
#########################################################

# Extract variables to be processed 
vars=c("Glucose","HDL_cholesterol","Total_cholesterol","whr","body_fat","bmi")

# Read in archived file that holds the order of IDs we want to achieve 
phen = read.table("bmi.phen",header=T)
ids=phen$FID

# Loop through each phenotype, reading it in and mean imputing before saving out 
for(i in vars){
  # Read in phenotype 
  phenos1=read.table(paste0(i,"_res_26092023.phen"),header=T) 
  # Add in the missing IDs 
  phenos1$IID = NULL # helps in merge step to prevent duplication of columns 
  phenos1=merge(phen[,c("FID","IID")],phenos1,by="FID",all.x=T)
  # Identify missingness
  names(phenos1)[3]=i # helps in automation
  message(paste0("There are ", length(which(is.na(phenos1[,i]))), " missing observations ", "for ", i))
  # Mean imputation step 
  phenos1[,i] <- meanimpute(phenos1[,i])
  # Match order of IDs with old files
  phenos1=phenos1[match(ids,phenos1$FID),]
  # Save out in correct format for BayesR+
  write.table(x = t(as.matrix(as.numeric(scale(phenos1[,i])))),file = paste0(i, ".csvphen"),quote = F, sep = ",", row.names = F, col.names = F)
  # Print to denote completion 
  message(paste0("Finished processing ", i))
} 

#########################################################
######## STEP 2 - PREPARATION OF COVARIATE FILES ########
#########################################################

# Read in WBCs
cov=read.csv("covariates.csv", header=T)
# Subset to columns of interest 
cov=cov[,c("Sample_Sentrix_ID","Bcell","CD8T","CD4T","NK","Mono","smokingScore")]
# Match order of IDs to other files 
cov=cov[match(ids,cov$Sample_Sentrix_ID),]
# Scale columns of interest 
cov[2:ncol(cov)]=apply(cov[,2:ncol(cov)],2,scale)
# Remove ID variables 
cov$Sample_Sentrix_ID=NULL
# Make dataframe without smokingScore for convenience in writing out step 
cov1=cov[,1:5]
# Save out file with just WBCs 
write.table(x = as.matrix(cov1),file = "wbcs_18413_covariates.csv" ,quote = F, sep = ",", row.names = F, col.names = F)
# Save out file with WBCs and smoking score 
write.table(x = as.matrix(cov),file = "wbcs_smk_18413_covariates.csv" ,quote = F, sep = ",", row.names = F, col.names = F)



#############################################
#### STEP 3 - BAYESR BASIC MODEL - WBCs #####
#############################################

for i in Phenotypes/Chain2/*.csvphen  

do 
A=$( echo $i | cut -d"/" -f3)
B=$( echo $A | cut -d"." -f1)

../../BayesRRcmd/src/brr --data-file Chain2/DNAm_18413_resid.csv --pheno $i --analysis-type preprocess --fixed_effects wbcs_18413_covariates_05062023.csv --fixedEffectNumber 5 --thread 12 --thread-spawned 12 --marker-cache --seed 1
../../BayesRRcmd/src/brr --data-file Chain2/DNAm_18413_resid.csv --pheno $i --fixed_effects wbcs_18413_covariates_05062023.csv --fixedEffectNumber 5 --analysis-type ppbayes --chain-length 10000 --burn-in 5000 --thin 5 --S "0.0001,0.001,0.01" --mcmc-samples Outputs_05062023/Chain2/${B}.csv --thread 12 --thread-spawned 12 --marker-cache --seed 1

echo $i
done 


############################################################
#### STEP 4 - BAYESR BASIC MODEL - WBCs + smokingScore #####
############################################################

for i in Phenotypes/*.csvphen  

do 
A=$( echo $i | cut -d"/" -f2)
B=$( echo $A | cut -d"." -f1)

../../BayesRRcmd/src/brr --data-file DNAm_18413_resid.csv --pheno $i --analysis-type preprocess --fixed_effects wbcs_smk_18413_covariates_05062023.csv --fixedEffectNumber 6 --thread 12 --thread-spawned 12 --marker-cache --seed 1
../../BayesRRcmd/src/brr --data-file DNAm_18413_resid.csv --pheno $i --fixed_effects wbcs_smk_18413_covariates_05062023.csv --fixedEffectNumber 6 --analysis-type ppbayes --chain-length 10000 --burn-in 5000 --thin 5 --S "0.0001,0.001,0.01" --mcmc-samples Outputs_05062023/${B}.csv --thread 12 --thread-spawned 12 --marker-cache --seed 1

echo $i
done 




