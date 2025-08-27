library(readr)

#read in episcores 
KORA_EpiScores_projected_in_LBC1936_WAVE1_17072023 <- read.csv("KORA_EpiScores_projected_in_LBC1936_WAVE1_17072023.csv")
EpiScores <- KORA_EpiScores_projected_in_LBC1936_WAVE1_17072023

#convert set, date and array to factor
EpiScores$set <- as.factor(EpiScores$set)
EpiScores$date <- as.factor(EpiScores$date)
EpiScores$array <- as.factor(EpiScores$array)


#check data is in a dataframe
class(EpiScores) # dataframe

#load lme4
library(lme4)

#make list of episcores
list_e <- colnames(EpiScores)[2:85]

#for loop that loops over each episcore adjusting for set, array and data
# takes the resdiuals and assigns them as the new episcores 
for(i in list_e){
  mod <- lmer(EpiScores[,i] ~ (1|set) + (1|array) + (1|date), 
            na.action = na.exclude, data = EpiScores)
  EpiScores[,i] <- resid(mod) 
}


# Rank-Inverse Based Normaliation 
library(bestNormalize)
for(i in colnames(EpiScores)[2:85]){ 
  EpiScores[,i]<- orderNorm(unlist(EpiScores[,i]))$x.t
}


#save out file
write.csv(EpiScores, file = "KORA_EpiScores_projected_Adjusted_in_LBC1936_WAVE1_17072023.csv")

