#read in episcores
episcore = readRDS("")

#read in test data
test = read.csv("testing_data_20.csv")


library(dplyr)

#join together
data = left_join(test, episcore, by = "Sample")

#correlate measure protein and episcore
cor(data$protein, data$EpiScore)

#test correlation measure protein and episcore
cor.test(data$protein, data$EpiScore)

