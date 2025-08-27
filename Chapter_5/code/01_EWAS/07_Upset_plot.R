library(dplyr)

bmi <- read.csv("bmi_EWAS_20methPCs_18411_significant.csv")

whr <- read.csv("whr_EWAS_20methPCs_18411_significant.csv")

body <- read.csv(#body_fat_EWAS_20methPCs_18411_significant.csv")

HDL <- read.csv("HDL_cholesterol_EWAS_20methPCs_18411_significant.csv")

Total <- read.csv("Total_cholesterol_20methPCs_18411_significant.csv")

Glucose <- read.csv("Glucose_EWAS_20methPCs_18411_significant.csv")

#select cpg and beta columns 
bmi <- bmi %>% select(Probe, b)
body <- body %>% select(Probe, b)
whr <- whr %>% select(Probe, b)
HDL <- HDL %>% select(Probe, b)
Total <- Total %>% select(Probe, b)
Glucose <- Glucose %>% select(Probe, b)

#rename phenotype columns
names(bmi)[2] <- "BMI (kg/m2)"
names(body)[2] <- "Body fat (%)"
names(whr)[2] <- "Waist-to-hip ratio"
names(HDL)[2] <- "HDL cholesterol (mmol/L)"
names(Total)[2] <- "Total cholesterol (mmol/L)"
names(Glucose)[2] <- "Glucose (mmol/L)"

#join results files
all_cpgs <- bmi %>% full_join(whr, by = "Probe")
all_cpgs <- all_cpgs %>% full_join(body, by = "Probe")
all_cpgs <- all_cpgs %>% full_join(HDL, by ="Probe")
all_cpgs <- all_cpgs %>% full_join(Total, by = "Probe")
all_cpgs <- all_cpgs %>% full_join(Glucose, by = "Probe")

#set CpG names ro rownames
row.names(all_cpgs) <- all_cpgs$Probe

#change values to 1 and "NA" to 0
all_cpgs[!is.na(all_cpgs)] <- 1
all_cpgs[is.na(all_cpgs)] <- 0

#set text options 
text_scale_options1 <- c(1, 1, 1, 1, 0.75, 1)
text_scale_options2 <- c(1.3, 1.3, 1, 1, 2, 0.75)
text_scale_options3 <- c(1.7, 1.5, 1.5, 1.5, 2, 1.5)

#setting colors
#this can also be done with hexadecimal
main_bar_col <- c("violetred4")
sets_bar_col <- c("turquoise4")
matrix_col <- c("slateblue4")
shade_col <- c("wheat4")


#set mb.ratio
mb_ratio1 <- c(0.55,0.45)


#set names are phenotype names
set_vars <- colnames(all_cpgs)[2:7]


#generate upset plot with UpSetR package
library(UpSetR)
setwd("")
tiff("Common_CpG_Upset_plot_20methPCs_18411.tiff", width = 16*300, height = 7*300, res=300)
      upset(all_cpgs, 
      sets = set_vars,
      mb.ratio = mb_ratio1, 
      mainbar.y.label = "Counts by Pattern of Conditions", 
      sets.x.label = "Counts by Condition",
      order.by = "freq",
      show.numbers = "yes",
      point.size = 2, 
      line.size = 1,
      text.scale=text_scale_options3,
      main.bar.color = main_bar_col,
      sets.bar.color = sets_bar_col,
      matrix.color = matrix_col,
      shade.color = shade_col,
      #nintersects= 28, 
      set_size.scale_max = 12400)
dev.off()
