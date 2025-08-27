library(readr)

#read in EpiScores
KORA_EpiScores_projected_in_LBC1921_WAVE1_24012023 <- read.csv("KORA_EpiScores_projected_in_LBC1921_WAVE1_24012023.csv")
EpiScores <- KORA_EpiScores_projected_in_LBC1921_WAVE1_24012023

#make set,date and arraye as factors
EpiScores$date <- as.factor(EpiScores$date)
EpiScores$array <- as.factor(EpiScores$array)

#check dataframe 
class(EpiScores) # dataframe

#load lme4 library 
library(lme4)

#make lsit of EpiScores
list_e <- colnames(EpiScores)[2:85]

#loop over EpiScores adjusting for array and date, extract residuals and assign to EpiScore
for(i in list_e){
  mod <- lmer(EpiScores[,i] ~ (1|array) + (1|date), 
              na.action = na.exclude, data = EpiScores)
  EpiScores[,i] <- resid(mod) 
}


# Rank-Inverse Based Normaliation 
library(bestNormalize)
for(i in colnames(EpiScores)[2:85]){ 
  EpiScores[,i]<- orderNorm(unlist(EpiScores[,i]))$x.t
}

#save out adjusted EpiScores
write.csv(EpiScores, file = "KORA_EpiScores_projected_Adjusted_in_LBC1921_WAVE1_29012023.csv")
