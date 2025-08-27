library(tidyverse)

# Dementia 

# Read in dementia 
dem <- read.csv("dementia_case_update_20241024.csv")

#check gp consent ids 
gp <- dem[which(dem$source == "gp"),]
dim(gp)
#[1] 88 15

#select id and data of diagnosis 
all <- dem[c(1,10)]

#rename columns
names(all) <- c("id", "first")

#prepare dementia file
all <- all[order(all$first),] # Check that the file is ordered by date (lowest to highest)
all <- all[which(all$first > 199001),] # Remove codes pre 199001
all <-all[-which(duplicated(all$id)),] # Remove any duplicate ID codes by taking the earliest date for that person 
dim(all) # 354
all <- all[which(all$first <= 202310), ] #subet to before censor date October 2023
dim(all)
#[1] 354  ## number of cases 

#save out dementia file
write.csv(all, "all_dementia_prepped.csv", row.names = F)

