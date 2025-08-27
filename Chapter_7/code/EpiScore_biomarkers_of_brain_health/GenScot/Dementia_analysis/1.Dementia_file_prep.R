library(tidyverse)

# Dementia 

# Read in dementia codesw
dem <- read.csv("Dementia_case_update_14Aug2023.csv")


# Summarise subgroup counts
table(dem$dementia) # 595
table(dem$ad) # 195
table(dem$vd) # 167
table(dem$dlb) # 1
table(dem$ftd) #7


# Check that all GP codes have given consent for GP sharing 
consent <- read.csv("2023-08-15_id_consent.csv")
consent <- as.data.frame(consent)
consent$Consent <- "YES"

dim(consent)
#[1] 7580    2

#join dementia codes and consent ids
dem <- left_join(dem, consent, by = "id")

#check gp consent ids 
gp <- dem[which(dem$source == "gp"),]
dim(gp)
#[1] 88 15

table(gp$Consent) 
#YES
#88

#subset to participants with dementia
all <- dem[which(dem$dementia == "1"),]

#select id and data of diagnosis 
all <- all[c(2,10)]

#rename columns
names(all) <- c("id", "first")

#prepare dementia file
all <- all[order(all$first),] # Check that the file is ordered by date (lowest to highest)
all <- all[which(all$first > 199001),] # Remove codes pre 199001
all <-all[-which(duplicated(all$id)),] # Remove any duplicate ID codes by taking the earliest date for that person (if ordered this should take first instance)
dim(all) # 334
all <- all[which(all$first <= 202204), ] #subet to before censor date April 2022
dim(all)
#[1] 334  ## number of cases 

#save out dementia file
write.csv(all, "all_dementia_prepped_15082023.csv", row.names = F)

