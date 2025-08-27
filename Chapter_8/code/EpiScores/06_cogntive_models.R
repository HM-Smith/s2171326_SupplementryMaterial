library(haven)

#read in cognitive data
df = read_sas("cogfac_1to10_10jul19.sas7bdat")

#rename column
names(df)[1] = "TWINNR"

#read in test data
test = read.csv("testing_data_20.csv")

#subset cognitive data to test data
df = df[df$TWINNR %in% test$TWINNR, ]

#select id and general cognitive tests
df = df[,c("TWINNR","tpcom5","tpcom6","tpcom8",
           "tpcom9","tpcom10")]

#put data into long format 
library(tidyr)
df_long <- pivot_longer(df,
                        cols = c("tpcom5","tpcom6","tpcom8",
                                 "tpcom9","tpcom10"),
                        names_to = "STUDY",
                        values_to = "G")

library(dplyr)

#rename testing appointments 
df_long$STUDY <- recode(df_long$STUDY,
                         "tpcom10" = "IPT10",
                         "tpcom5"  = "IPT5",
                         "tpcom6"  = "IPT6",
                         "tpcom8"  = "IPT8",
                         "tpcom9"  = "IPT9")

#subset to complete daya
df_long = df_long[complete.cases(df_long),]

#check unique participants
length(unique(df_long$TWINNR))

#merge datasets
cog = left_join(df_long, test, by = c("TWINNR", "STUDY"))

#read in episcore
episcore = readRDS("")

#join cogntive data and episcore
cog = left_join(cog, episcore, by = c("Sample"))

#subset to complete data
cog = cog[complete.cases(cog), ]



library(lme4)
library(lmerTest)

#run mixed effects models for episcore and G
mod <- lmer(scale(G) ~ scale(episcore) + 
              scale(AGE) + 
              as.factor(SEX) +
              (1|PAIRID) +
              (1|TWINNR), data = cog)

#summary episcore model
summary(mod)

#run mixed effects models for protein and G
mod2 <- lmer(scale(G) ~ scale(protein) + 
              scale(AGE) + 
              as.factor(SEX) +
              (1|PAIRID) +
              (1|TWINNR), data = cog)

#summary protein model
summary(mod2)