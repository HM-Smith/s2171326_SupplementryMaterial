library(dplyr)
library(lme4)

#read in episcores
EpiScores <- read.csv("LBC1936_with_all_lipid_EpiScores.csv")

#set date and array as factors
EpiScores$date <- as.factor(EpiScores$date)
EpiScores$array <- as.factor(EpiScores$array)


#make list of episcores
list_e <- colnames(EpiScores)[109:114]

#for loop that loops over each episcore adjusting for array and data
# takes the resdiuals and assigns them as the new episcores 
for(i in list_e){
  mod <- lmer(EpiScores[,i] ~ (1|array) + (1|date), 
              na.action = na.exclude, data = EpiScores)
  EpiScores[,i] <- resid(mod) 
}

#save out adjusted episcores
write.csv(EpiScores, file = "LBC1936_with_all_lipid_EpiScores_array_date_adj.csv")
