#######################
### read in IQ data ###
#######################

#read in cogntive data 
cog = read.csv("cognitive.csv", header=T)

#add together cognitive test scores to create logical memory variable 
cog$LM <- cog$logical_mem_1 + cog$logical_mem_2

##############################################
### recode 0 to NA for each cognitive test ###
##############################################
table(cog$LM==0)
#FALSE  TRUE
#21239     2

table(cog$digit_symbol==0)
#FALSE  TRUE
#21259     7

table(cog$vocabulary==0)
#FALSE  TRUE
#21115     7

table(cog$verbal_total==0)
#FALSE  TRUE
#21250     4

table(is.na(cog))
#FALSE   TRUE
#467747   5781


cog$LM[cog$LM==0] <- NA
cog$digit_symbol[cog$digit_symbol==0] <- NA
cog$vocabulary[cog$vocabulary==0] <- NA
cog$verbal_total[cog$verbal_total==0] <- NA

table(is.na(cog))
#FALSE   TRUE
#467727   5801

5801 - 5781 #20 - as expected from no. of zeros in cog tests 




############################################################################
### recode points +- 3.5 SDs from the mean to NA for each cognitive test ### 
############################################################################

low = mean(cog$LM, na.rm=T) - 3.5*sd(cog$LM, na.rm=T)
high = mean(cog$LM, na.rm=T) + 3.5*sd(cog$LM, na.rm=T)
table(cog$LM < low | cog$LM > high)
#FALSE  TRUE
#21221    18
cog$LM[cog$LM < low | cog$LM > high] <- NA

low = mean(cog$verbal_total, na.rm=T) - 3.5*sd(cog$verbal_total, na.rm=T)
high = mean(cog$verbal_total, na.rm=T) + 3.5*sd(cog$verbal_total, na.rm=T)
table(cog$verbal_total < low | cog$verbal_total > high)
#FALSE  TRUE
#21214    36
cog$verbal_total[cog$verbal_total < low | cog$verbal_total > high] <- NA

low = mean(cog$digit_symbol, na.rm=T) - 3.5*sd(cog$digit_symbol, na.rm=T)
high = mean(cog$digit_symbol, na.rm=T) + 3.5*sd(cog$digit_symbol, na.rm=T)
table(cog$digit_symbol < low | cog$digit_symbol > high)
#FALSE  TRUE
#21234    25
cog$digit_symbol[cog$digit_symbol < low | cog$digit_symbol > high] <- NA

low = mean(cog$vocabulary, na.rm=T) - 3.5*sd(cog$vocabulary, na.rm=T)
high = mean(cog$vocabulary, na.rm=T) + 3.5*sd(cog$vocabulary, na.rm=T)
table(cog$vocabulary < low | cog$vocabulary > high)
#FALSE  TRUE
#21052    63
cog$vocabulary[cog$vocabulary < low | cog$vocabulary > high] <- NA


18 + 36 + 25 + 63 # 142
table(is.na(cog))
#FALSE   TRUE
#467585   5943

5943 - 5801 # 142, as expected from the no. of outliers in cog tests


# rename id column 
colnames(cog)[1] <- "Sample_Name"

#read in episcores
GS_Episcores <- read.csv("KORA_Episcores_projected_adjusted_in_GS20K_19012023.csv")

#join cognitive data and episcores
GS_cog <- cog %>% left_join(GS_Episcores, by = "Sample_Name")

#recode sex from characters to binary 
GS_cog <- mutate(GS_cog, sex=recode(sex, 'M' = '1', 'F' = '2'))

#set sex as factor
GS_cog$sex <- as.factor(GS_cog$sex)

#save out cognitive data with episcores
write.csv(GS_cog, file = "GS20K_Cognitive_KORA_adjust_EpiScores_combined_20012023.csv")

