###########################################################################################################

### GS EpiScore prep

###########################################################################################################

library(tidyverse)

# Load EpiScores in 20k projections and combine
W1W3 <- read.csv("Projected_Episcores/EpiScore_projections_GS_9537.csv", check.names = F)
W4 <- read.csv("Projected_Episcores/EpiScore_projections_W4_8877_220221.csv", check.names = F)
names(W4)[1] <- "ID"
colnames(W1W3) == colnames(W4)
scores <- rbind(W1W3, W4) 
names(scores)[1] <- "Sample_Sentrix_ID"
scores <- scores %>% select(1, 27:ncol(scores))

#load sample ID file 
target <- readRDS("GS20k_Targets.rds")

# Add sample ID information for joining 
scores <- left_join(scores, target, by = "Sample_Sentrix_ID") 

class(scores) #dataframe

#loop over episcores and adjust scores for set and batch then extract residuals 
library(lme4)
list_e <- colnames(scores)[2:85]
for(i in list_e){
  mod <- lmer(scores[,i] ~ (1|Set) + (1|Batch), 
              na.action = na.exclude, data = scores)
  scores[,i] <- resid(mod) 
}


# Transform episcores  
library(bestNormalize)
for(i in colnames(scores)[2:85]){ 
  scores[,i]<- orderNorm(scores[,i])$x.t # Rank-Inverse Based Normaliation
}

# Save episcores prepped 
write.csv(scores, file = "KORA_Episcores_projected_adjusted_in_GS20K_19012023.csv", row.names = F)
