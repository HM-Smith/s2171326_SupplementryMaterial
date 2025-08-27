library(dplyr)
library(readxl)
library(ggplot2)
library(patchwork)

#read in cogntive function assocs full model 
data <- read.csv("G_intercept_No_Domain_LBC1936_assocs_lipid_EpiScores_with_ids_full_no_yrs_edu.csv")

#print traits
data$SeqId
#[1] "bld_choles_w1"  "bld_hdlchol_w1" "BMI_Epi"        "bmi_w1"
#[5] "Body_fat_Epi"   "Glucose_Epi"    "HDL_Epi"        "Total_chol_Epi"
#[9] "WHR_Epi"


#rename traits
data$Trait <- c("Total cholesterol", "HDL cholesterol", "BMI", "BMI", "Body fat", "Glucose", "HDL cholesterol", "Total cholesterol", "WHR")

#create "measured or episcore" column 
data$Type <- c("Measured","Measured", "EpiScore", "Measured", "EpiScore","EpiScore","EpiScore", "EpiScore", "EpiScore")

#order by effect size for plotting
data$SeqId = factor(data$SeqId, levels=unique(data$SeqId[order(data$beta)]))

#make a column that indicates level of significant or non-significance
data$Significance <- ifelse(data$FDR < 0.05, "FDR significant", ifelse(data$P < 0.05, "Nominally significant", "Non significant"))


#plot results as forrest plot and save as p1
p1 <- ggplot(data, aes(y = beta, x = Trait, color = Type, shape = Significance)) +
  geom_point(size = 2.5, position=position_dodge(0.5)) +
  geom_errorbar(aes(ymax = ci.upper, ymin = ci.lower), size = .6, width =
                  0, position=position_dodge(0.5)) +
  geom_hline(aes(yintercept = 0), size = .25, linetype = "dashed") +
  #coord_trans(y = scales:::exp_trans(10)) + 
  theme_classic() +
  theme(legend.title=element_blank())+
  #theme(panel.grid.minor = element_blank()) +
  ylab("Standardised Beta [95% CI]") +
  xlab("Trait") +
  #ggtitle("General cognitive function in LBC1936") +
  theme(plot.title = element_text(hjust = 0.5, size = 12))+
  #guides(color = "none", size = "none")+
  #guides(shape = guide_legend(override.aes = list(size=5))) + 
  coord_flip()

#read in variance explain in cog function results
data <- read_excel("Cognitive_R2_measure_episcore.xlsx")

#select only required columns
data <- data[-1,]

#plot var explain results as bar chart and save as p2
p2 <- ggplot(data, aes(x=Trait, y=`Incremental R2`, fill = Group)) +
  geom_bar(stat="identity",  position=position_dodge()) +
  ylab(expression("Incremental R"^2)) +
  theme_classic()+
  theme(legend.title=element_blank())


#combine plot 1 and 2 and multi panel plot with labels (A and B)
setwd("")
tiff("figure2_multi_panel.tiff", width =20*300, height = 18*300, res = 300 )

(p1 + p2) + 
  plot_annotation(tag_levels = 'A') & 
  theme(text = element_text(family = "Arial"))

dev.off()

#################################################################################################
####################### plot of basic cog model for supplementary ################################

library(dplyr)

#read in basic results file
data <- read.csv("G_intercept_No_Domain_LBC1936_assocs_lipid_EpiScores_with_ids_basic.csv")

#print traits
data$SeqId
#[1] "bld_choles_w1"  "bld_hdlchol_w1" "BMI_Epi"        "bmi_w1"
#[5] "Body_fat_Epi"   "Glucose_Epi"    "HDL_Epi"        "Total_chol_Epi"
#[9] "WHR_Epi"


#rename trait columns
data$Trait <- c("Total cholesterol","HDL cholesterol", "BMI", "BMI", "Body fat", "Glucose", "HDL cholesterol","Total cholesterol", "WHR")

#make "measured or episcore" column
data$Type <- c("Measured","Measured", "EpiScore", "Measured", "EpiScore","EpiScore","EpiScore","EpiScore", "EpiScore")

#order by effect size for plotting
data$SeqId = factor(data$SeqId, levels=unique(data$SeqId[order(data$beta)]))


library(ggplot2)

#plot results as forrest plot
setwd("")
tiff("Cognitive_function_basic_models.tiff", width =5.5*300, height = 6*300, res = 300)
ggplot(data, aes(y = beta, x = Trait, color = Type)) +
  geom_point(size = 2.5, position=position_dodge(0.5)) +
  geom_errorbar(aes(ymax = ci.upper, ymin = ci.lower), size = .6, width =
                  0, position=position_dodge(0.5)) +
  geom_hline(aes(yintercept = 0), size = .25, linetype = "dashed") +
  #coord_trans(y = scales:::exp_trans(10)) + 
  theme_classic() +
  theme(legend.title=element_blank())+
  #theme(panel.grid.minor = element_blank()) +
  ylab("Standardised Beta [95% CI]") +
  xlab("Trait") +
  #ggtitle("General cognitive function in LBC1936") +
  theme(plot.title = element_text(hjust = 0.5, size = 12))+
  #guides(color = "none", size = "none")+
  #guides(shape = guide_legend(override.aes = list(size=5))) + 
  coord_flip()
dev.off()

