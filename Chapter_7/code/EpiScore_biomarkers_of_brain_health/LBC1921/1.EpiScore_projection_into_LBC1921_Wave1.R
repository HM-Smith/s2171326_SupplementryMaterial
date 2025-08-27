############################################################################################

## Projections in LBC1936
library(tidyverse)

# Read in LBC methylation file, methylation is uncorrected for any technical
# covariates and is in Beta value format
load("Beta_3525_norm_bgcorrect_0.001BetaThreshold_probefilter.RObject") 

# Read in EPIC annotation file
anno <- readRDS("EPIC_AnnotationObject_df.rds")

# Subset annotation file to probes common to 450k and EPIC array
common_anno <- anno[which(anno$Methyl450_Loci == "TRUE"),]

# Subset methylation data file to those in the common annotation file 
dat1 <- dat[rownames(dat) %in% rownames(common_anno),] 

# > dim(dat1)
# [1] 428489   3525

# Read in target file with matched IDs for methylation data 
target <- read.csv("target_QC_age_sex_date.csv")

# Filter to wave 1
dat3 <- target %>% filter(WAVE == "1") # 1342 people

# Filter to LBC36
dat3  <- dat3 %>% filter(cohort == "LBC21") # 906 people 

# Subset DNAm data to these individuals from W1 of LBC36 only 
dat1 <- dat1[,which(colnames(dat1) %in% dat3$Basename)]

dim(dat1)
#[1] 428489    436


# Transpose data so CpGs are columns 
dat1 <- t(dat1)

# Replace NA values in methylation data with imputed means (means for each CpG column imputed)
for(t in 1:ncol(dat1)){
  dat1[is.na(dat1[,t]), t] <- mean(dat1[,t], na.rm = TRUE)
}

# Assign variable for data to the LBC DNAm dataset 
data = dat1

# Transpose back so CpGs are rows again 
data <- t(data)



## Process weights files 
cpgs_KORA <- read.csv("KORA_weights_finalised_120121.csv")
cpgs_KORA$Mean_Beta_Value <- 0
cpgs_KORA <- cpgs_KORA[c(1,2,14,7)] # dim(cpgs_KORA) 8374, 4
names(cpgs_KORA) <- c("CpG_Site", "Coefficient", "Mean_Beta_Value", "Predictor")


cpgs <- cpgs_KORA

#check unique CpGs
un <- unique(cpgs$Predictor)
un <- as.data.frame(un) # 84 episcores
names(un) <- "Protein"



## Check if Data needs to be Transposed
if(ncol(data) > nrow(data)){
  message("It seems that individuals are rows - data will be transposed!")
  data<-t(data)}

## Subset CpG sites to those present on list for predictors 
coef=data[intersect(rownames(data), cpgs$CpG_Site),] # 7407 cpgs of 8374 present (KORA only)


## Check if Beta or M Values, convert to beta values if required
m_to_beta <- function (val) 
{
  beta <- 2^val/(2^val + 1)
  return(beta)}

coef<-if((range(coef,na.rm=T)> 1)[[2]] == "TRUE") { message("Suspect that M Values are present. Converting to Beta Values");m_to_beta(coef) } else { message("Suspect that Beta Values are present");coef}

## Scale Data if Needed 
ids = colnames(coef)
scaled <- apply(coef, 1, function(x) sd(x,na.rm = T)) 

coef <-  if(range(scaled)[1] == 1 & range(scaled)[2] == 1) { 
  coef
} else { 
  coef_scale <- apply(coef, 1, scale)
  coef_scale <- t(coef_scale)
  coef_scale <- as.data.frame(coef_scale)
  colnames(coef_scale) <- ids
  coef_scale
} 


## Identify CpGs missing from input dataframe, include them and provide values as mean methylation value at that site
## 170 sites are missing for LBC1936 WAVE 1
coef <- if(nrow(coef) == 10270) { message("All sites present"); coef } else if(nrow(coef)==0){ 
  message("There Are No Necessary CpGs in The dataset - All Individuals Would Have Same Values For Predictors. Analysis Is Not Informative!")
} else { 
  missing_cpgs = cpgs[-which(cpgs$CpG_Site %in% rownames(coef)),c("CpG_Site","Mean_Beta_Value")]
  message(paste(length(unique(missing_cpgs$CpG_Site)), "unique sites are missing - add to dataset with mean Beta Value from Training Sample", sep = " "))
  mat = matrix(nrow=length(unique(missing_cpgs$CpG_Site)),ncol = ncol(coef))
  row.names(mat) <- unique(missing_cpgs$CpG_Site)
  colnames(mat) <- colnames(coef) 
  mat[is.na(mat)] <- 1
  missing_cpgs1 <- if(length(which(duplicated(missing_cpgs$CpG_Site))) > 1) { 
    missing_cpgs[-which(duplicated(missing_cpgs$CpG_Site)),]
  } else {missing_cpgs
  }  
  ids = unique(row.names(mat))
  missing_cpgs1 = missing_cpgs1[match(ids,missing_cpgs1$CpG_Site),]
  mat=mat*missing_cpgs1$Mean_Beta_Value
  coef=rbind(coef,mat) } 


## Convert NAs to Mean Value for all individuals across each probe 
na_to_mean <-function(methyl){
  methyl[is.na(methyl)]<-mean(methyl,na.rm=T)
  return(methyl)
}

coef <- t(apply(coef,1,function(x) na_to_mean(x)))


#loop to calculate EpiScores 
loop = unique(cpgs$Predictor)
out <- data.frame()
for(i in loop){ 
  tmp=coef[intersect(row.names(coef),cpgs[cpgs$Predictor %in% i,"CpG_Site"]),]
  tmp_coef = cpgs[cpgs$Predictor %in% i, ]
  if(nrow(tmp_coef) > 1) { 
    tmp_coef = tmp_coef[match(row.names(tmp),tmp_coef$CpG_Site),]
    out[colnames(coef),i]=colSums(tmp_coef$Coefficient*tmp)
  } else {
    tmp2 = as.matrix(tmp)*tmp_coef$Coefficient 
    out[colnames(coef),i] = tmp2[,1]
  }
}

#make rownames Basename column 
out$Basename <- row.names(out)

# merge target in with episcores 
out_final <- out %>% left_join(target, by="Basename")

dim(out_final)
#[1] 436 102


## Save file 
write.csv(out_final, "KORA_EpiScores_projected_in_LBC1921_WAVE1_24012023.csv")