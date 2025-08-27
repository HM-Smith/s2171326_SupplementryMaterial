library(tidyverse)
library(ggpattern)
library(readxl)

#read in elnet episcores variance explained table 
HELIOS_LBC1936_combined_incremental_Rsquared_results <- read_excel("HELIOS_LBC1936_Elnet_combined_incremental_Rsquared_results.xlsx")
data <- HELIOS_LBC1936_combined_incremental_Rsquared_results

#set order of Groups so LBC1936 and HELIOS-all appear next to eachother in graph
data$Group = factor(data$Group, levels = c('LBC1936', 'HELIOS - all', 'HELIOS - Chinese', 'HELIOS - Indian', 'HELIOS - Malay'))

#plot variance explained results as bar chart 
tiff("IncrementalR2_plot_Elnet.tiff", width = 14*300, height = 9*300, res = 300)
ggplot(data, aes(x=`Trait`, y=`Incremental R2`, fill = Group)) +
  geom_bar(stat = "identity", position = position_dodge(), color = "black")+
  ylab(expression("Incremental R"^2)) +
  xlab("Trait")+
  scale_x_discrete(labels = c("BMI" = expression("BMI (kg/m"^2*")"), "Body fat %" = "Body fat (%)", "HDL cholesterol" = "HDL cholesterol (mmol/L)",
                              "Total cholesterol" = "Total cholesterol (mmol/L)", "WHR" = "WHR")) +
  theme_classic() + 
  theme(text=element_text(size=15)) +
  ylim(min = 0, max = 30)

dev.off()

#remove dataset to re-run with BayesR episcore results
rm(HELIOS_LBC1936_combined_incremental_Rsquared_results)
rm(data)

###########################################################################################
##########################################################################################

#read in BayesR episcore variance explained results
HELIOS_LBC1936_combined_incremental_Rsquared_results <- read_excel("HELIOS_LBC1936_BayesR_combined_incremental_Rsquared_results.xlsx")
data <- HELIOS_LBC1936_combined_incremental_Rsquared_results

#set order of Groups so LBC1936 and HELIOS-all appear next to eachother in graph
data$Group = factor(data$Group, levels = c('LBC1936', 'HELIOS - all', 'HELIOS - Chinese', 'HELIOS - Indian', 'HELIOS - Malay'))

#plot variance explained results as a bar chart
tiff("IncrementalR2_plot_BayesR.tiff", width = 14*300, height = 9*300, res = 300)
ggplot(data, aes(x=`Trait`, y=`Incremental R2`, fill = Group)) +
  geom_bar(stat = "identity", position = position_dodge(), color = "black")+
  ylab(expression("Incremental R"^2)) +
  xlab("Trait")+
  scale_x_discrete(labels = c("BMI" = expression("BMI (kg/m"^2*")"), "Body fat %" = "Body fat (%)", "HDL cholesterol" = "HDL cholesterol (mmol/L)",
                              "Total cholesterol" = "Total cholesterol (mmol/L)", "WHR" = "WHR")) +
  theme_classic() + 
  theme(text=element_text(size=15)) +
  ylim(min = 0, max = 30)
dev.off()
