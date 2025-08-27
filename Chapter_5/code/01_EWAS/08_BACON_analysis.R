library("bacon")
library("dplyr")

#read in results 
BMI <- read.csv("bmi_EWAS_20methPCs_18411_filtered_CpGs.csv")

#make bacon object
BMI_INF <- bacon(
  effectsizes = BMI$b,
  standarderrors = BMI$se)

#extract corrected betas 
BMI$beta_INF <- es(BMI_INF)

#extract corrected p values 
BMI$P_INF <- pval(BMI_INF)

#extract corrected standard errors 
BMI$se_INF <- se(BMI_INF)

#save out corrected results 
write.csv(BMI, file = "bmi_EWAS_20methPCs_18411_filtered_CpGs_INF_Adj_08082024.csv")

#filter to significant hits
bmi_sig <- BMI[BMI$P_INF < 3.6e-8,]

#check num of sig hits
dim(bmi_sig)
#[1] 3710   13

#save out sig hits 
write.csv(bmi_sig, file = "bmi_EWAS_20methPCs_18411_significant_INF_Adj.csv")



#read in results 
body <- read.csv("body_fat_EWAS_20methPCs__18411_filtered_CpGs.csv")

#make bacon object
body_INF <- bacon(
  effectsizes = body$b,
  standarderrors = body$se)

#extract corrected betas 
body$beta_INF <- es(body_INF)
#extract corrected p values 
body$P_INF <- pval(body_INF)
#extract corrected standard errors 
body$se_INF <- se(body_INF)

#save out corrected results 
write.csv(body, file = "body_fat_EWAS_20methPCs__18411_filtered_CpGs_INF_Adj.csv")

#filter to significant hits
body_sig <- body[body$P_INF < 3.6e-8,]

#check num of sig hits
dim(body_sig)
#[1] 4003   13

#save out sig hits 
write.csv(body_sig, file = "body_fat_EWAS_20methPCs_18411_significantINF_Adj_080.csv")


#read in results 
glucose <- read.csv("Glucose_EWAS_20methPCs_18411_filtered_CpGs.csv")

#make bacon object
glucose_INF <- bacon(
  effectsizes = glucose$b,
  standarderrors = glucose$se)

#extract corrected betas 
glucose$beta_INF <- es(glucose_INF)
#extract corrected p values 
glucose$P_INF <- pval(glucose_INF)
#extract corrected standard errors 
glucose$se_INF <- se(glucose_INF)

#save out corrected results 
write.csv(glucose, file = "Glucose_EWAS_20methPCs_18411_filtered_CpGs_INF_Adj.csv")

#filter to significant hits
glucose_sig <- glucose[glucose$P_INF < 3.6e-8,]

#check num of sig hits
dim(glucose_sig)
#[1] 206  13

#save out sig hits 
write.csv(glucose_sig, file = "Glucose_EWAS_20methPCs_18411_significant_INF_Adj.csv" )


#read in results 
HDL <- read.csv("HDL_cholesterol_EWAS_20methPCs_18411_filtered_CpGs.csv")

#make bacon object
HDL_INF <- bacon(
  effectsizes = HDL$b,
  standarderrors = HDL$se)

#extract corrected betas 
HDL$beta_INF <- es(HDL_INF)
#extract corrected p values 
HDL$P_INF <- pval(HDL_INF)
#extract corrected standard errors 
HDL$se_INF <- se(HDL_INF)

#save out corrected results 
write.csv(HDL, file = "HDL_cholesterol_EWAS_20methPCs_18411_filtered_CpGs_INF_Adj.csv")

#filter to significant hits
HDL_sig <- HDL[HDL$P_INF < 3.6e-8,]

#check num of sig hits
dim(HDL_sig)
#[1] 4390   13

#save out sig hits 
write.csv(HDL_sig, file = "HDL_cholesterol_EWAS_20methPCs_18411_significant_INF_Adj.csv")


#read in results 
total <- read.csv("Total_cholesterol_EWAS_20methPCs_18411_filtered_CpGs.csv")

#make bacon object
total_INF <- bacon(
  effectsizes = total$b,
  standarderrors = total$se)

#extract corrected betas 
total$beta_INF <- es(total_INF)
#extract corrected p values 
total$P_INF <- pval(total_INF)
#extract corrected standard errors 
total$se_INF <- se(total_INF)

#save out corrected results 
write.csv(total, file = "Total_cholesterol_EWAS_20methPCs_18411_filtered_CpGs_INF_Adj.csv")

#filter to significant hits
total_sig <- total[total$P_INF < 3.6e-8,]

#check num of sig hits
dim(total_sig)
#[1] 877  13

#save out sig hits 
write.csv(total_sig, file = "Total_cholesterol_20methPCs_18411_significant_INF_Adj.csv")

#read in results 
WHR <- read.csv("whr_EWAS_20methPCs_18411_filtered_CpGs.csv")

#make bacon object
WHR_INF <- bacon(
  effectsizes = WHR$b,
  standarderrors = WHR$se)


#extract corrected betas 
WHR$beta_INF <- es(WHR_INF)
#extract corrected p values 
WHR$P_INF <- pval(WHR_INF)
#extract corrected standard errors 
WHR$se_INF <- se(WHR_INF)

#save out corrected results 
write.csv(WHR, file = "whr_EWAS_20methPCs_18411_filtered_CpGs_INF_Adj.csv")

#filter to significant hits
WHR_sig <- WHR[WHR$P_INF < 3.6e-8,]

#check num of sig hits
dim(WHR_sig)
#[1] 1421   13

#save out sig hits 
write.csv(WHR_sig, file = "whr_EWAS_20methPCs_18411_significant_INF_Adj.csv")
