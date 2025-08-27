library(dplyr)

#read in protein EWAS results 
df = read.csv("")

#filter to sig threshold
df_sig <- df %>% filter(P.value < 3.6e-8)

#filter to nominal significance
df_nom <- df %>% filter(P.value < 0.05)


write.csv(df_sig, file = "")
write.csv(df_nom, file = "")

