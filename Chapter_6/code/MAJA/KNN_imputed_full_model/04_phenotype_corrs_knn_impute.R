data = read.csv("complete_data_KNN_imputed.csv")

#select contious varaibles
cont <- data[,c(3:8,448:451, 454, 457)]

names(cont)

# [1] "digit_symbol"             "verbal_total"
# [3] "vocabulary"               "LM"
# [5] "alzheimers_scorefile_SUM" "cognition_scorefile_SUM"
# [7] "rank"                     "bmi"
# [9] "units"                    "age"
# [11] "years"                  "smokingScore"

#rename columns 
names(cont) <-c("digit symbol", "verbal fluency","vocabulary",
                "logical memory","AD PRS", "cognitive PRS", 
                "SMID", "BMI",
                "alcohol intake", "age",
                "years","smoking score")

#create correlation matrix
res2 <- cor(cont, use="pairwise.complete.obs")

#make correlation plot
library(corrplot)
setwd("")
tiff("outcomes_covar_corr_KNN_IMPUTED.tiff", width = 8*300, height = 9*300, res = 300)

corrplot(res2, type = "lower", order = "hclust", tl.col = "black", addCoef.col = 'black', method = "color", number.digits = 2)
dev.off()