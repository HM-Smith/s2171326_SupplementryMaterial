library(tidyverse)
library(ggplot2)

GenScot_FDR <- read.csv("/Cluster_Filespace/Marioni_Group/Hannah/Biomarkers_brain_health_LBC_GS/Results/GenScot/Dementia/Result_analysis/FDR_significant/Coxme_Dementia_GenScot_84_episcores_FDR_sig_basic_16082023.csv")

dim(GenScot_FDR)
#[1] 13 18

LBC1936_Nom_sig <- read.csv("/Cluster_Filespace/Marioni_Group/Hannah/Biomarkers_brain_health_LBC_GS/Results/LBC1936/Dementia/Result_analysis/nominally_significant/coxPH_Dementia_LBC1936_84_EpiScores_basic_model_nominally_sig_20072023.csv")

LBC1936_Nom_sig_common <- LBC1936_Nom_sig[LBC1936_Nom_sig$SeqId %in% GenScot_FDR$SeqId, ]

dim(LBC1936_Nom_sig_common)
#[1]  0 18


LBC1936_all <- read.csv("/Cluster_Filespace/Marioni_Group/Hannah/Biomarkers_brain_health_LBC_GS/Results/LBC1936/Dementia/Result_analysis/results_with_ids/coxPH_Dementia_LBC1936_84_EpiScores_basic_model_WITH_IDS_20072023.csv")
LBC1936_common <- LBC1936_all[LBC1936_all$SeqId %in% GenScot_FDR$SeqId,]

dim(LBC1936_common)
#[1]  13 18

LBC1936_common$Cohort <- rep("LBC1936", 13)
LBC1936_common$Model <- rep("coxPH", 13)
LBC1936_common$Significance <- rep("Non significant", 13)


GenScot_FDR$Cohort <- rep("GenScot", 13)
GenScot_FDR$Model <- rep("coxme", 13)
GenScot_FDR$Significance <- rep("FDR significant", 13)


GenScot_comprisk_Nom_sig <- read.csv("/Cluster_Filespace/Marioni_Group/Hannah/Biomarkers_brain_health_LBC_GS/Results/GenScot/Dementia/Result_analysis/nominally_significant/CompRisk_Dementia_GenScot_fineandgray_84_episcores_nominally_sig_basic_16082023.csv")

GenScot_comprisk_Nom_sig_common <- GenScot_comprisk_Nom_sig[GenScot_comprisk_Nom_sig$SeqId %in% GenScot_FDR$SeqId, ]

dim(GenScot_comprisk_Nom_sig_common)
#[1]  3 15

names(GenScot_comprisk_Nom_sig_common)[5] <- "Hazard_Ratio"
GenScot_comprisk_Nom_sig_common$Cohort <- rep("GenScot", 3)
GenScot_comprisk_Nom_sig_common$Model <- rep("Competing Risk", 3)
GenScot_comprisk_Nom_sig_common$Significance <- rep("Nominally significant", 3)


GenScot_all_comprisk <- read.csv("/Cluster_Filespace/Marioni_Group/Hannah/Biomarkers_brain_health_LBC_GS/Results/GenScot/Dementia/Result_analysis/results_with_ids/CompRisk_Dementia_GenScot_fineandgray_84_episcores_with_ids_basic_16082023.csv")

GenScot_all_comprisk_common <- GenScot_all_comprisk[GenScot_all_comprisk$SeqId %in% GenScot_FDR$SeqId & !GenScot_all_comprisk$SeqId %in% GenScot_comprisk_Nom_sig_common$SeqId,]

dim(GenScot_all_comprisk_common)
#[1]  10 15

names(GenScot_all_comprisk_common)[5] <- "Hazard_Ratio"
GenScot_all_comprisk_common$Cohort <- rep("GenScot", 10)
GenScot_all_comprisk_common$Model <- rep("Competing Risk", 10)
GenScot_all_comprisk_common$Significance <- rep("Non significant", 10)


LBC1936_comprisk_Nom_sig <- read.csv("/Cluster_Filespace/Marioni_Group/Hannah/Biomarkers_brain_health_LBC_GS/Results/LBC1936/Dementia/Result_analysis/nominally_significant/compRisk_Dementia_LBC1936_84_EpiScores_nominally_sig_basic_model_20072023.csv")

LBC1936_comprisk_Nom_sig_common <- LBC1936_comprisk_Nom_sig[LBC1936_comprisk_Nom_sig$SeqId %in% GenScot_FDR$SeqId,]

dim(LBC1936_comprisk_Nom_sig_common)
#[1]  0 15

LBC1936_comprisk_all <- read.csv("/Cluster_Filespace/Marioni_Group/Hannah/Biomarkers_brain_health_LBC_GS/Results/LBC1936/Dementia/Result_analysis/results_with_ids/compRisk_Dementia_LBC1936_84_EpiScores_with_ids_basic_model_20072023.csv")

LBC1936_comprisk_all_common <- LBC1936_comprisk_all[LBC1936_comprisk_all$SeqId %in% GenScot_FDR$SeqId,]

dim(LBC1936_comprisk_all_common)
#[1]  13 15

names(LBC1936_comprisk_all_common)[5] <- "Hazard_Ratio"
LBC1936_comprisk_all_common$Cohort <- rep("LBC1936", 13)
LBC1936_comprisk_all_common$Model <- rep("Competing Risk", 13)
LBC1936_comprisk_all_common$Significance <- rep("Non significant", 13)


GenScot_FDR <- GenScot_FDR %>% select(SeqId, Hazard_Ratio, ci.lower, ci.upper, Gene_Name, Cohort, Significance, Model)
LBC1936_common <- LBC1936_common %>% select(SeqId, Hazard_Ratio, ci.lower, ci.upper, Gene_Name, Cohort, Significance, Model)
GenScot_comprisk_Nom_sig_common <- GenScot_comprisk_Nom_sig_common %>% select(SeqId, Hazard_Ratio, ci.lower, ci.upper, Gene_Name, Cohort, Significance, Model)
GenScot_all_comprisk_common <- GenScot_all_comprisk_common %>% select(SeqId, Hazard_Ratio, ci.lower, ci.upper, Gene_Name, Cohort, Significance, Model)
LBC1936_comprisk_all_common <- LBC1936_comprisk_all_common %>% select(SeqId, Hazard_Ratio, ci.lower, ci.upper, Gene_Name, Cohort, Significance, Model)

all<- rbind(GenScot_FDR, LBC1936_common)
all <- rbind(all, GenScot_comprisk_Nom_sig_common)
all <- rbind(all, GenScot_all_comprisk_common)
all <- rbind(all, LBC1936_comprisk_all_common)

all$Gene_Name = factor(all$Gene_Name, levels=unique(all$Gene_Name[order(all$Hazard_Ratio)]))

setwd("/Cluster_Filespace/Marioni_Group/Hannah/Biomarkers_brain_health_LBC_GS/Results/combined_results/Dementia/")
pdf("Dementia_plot_supp_06112023.pdf")
ggplot(all, aes(y = Hazard_Ratio, x = Gene_Name, color = Cohort, shape = Significance, linetype = Model)) +
  geom_point(size = 2.5, position=position_dodge(0.5)) +
  geom_errorbar(aes(ymax = ci.upper, ymin = ci.lower),size = .5, width =
                  0, position=position_dodge(0.5)) +
  scale_color_manual(values = c('#F8766D', '#619CFF')) +
  scale_shape_manual(values = c(17,16,0)) +
  geom_hline(aes(yintercept = 1), size = .25, linetype = "dashed") +
  #coord_trans(y = scales:::exp_trans(10)) + 
  theme_classic() +
  #theme(panel.grid.minor = element_blank()) +
  ylab("Hazard Ratio [95% CI]") +
  xlab("EpiScore") +
  ggtitle("Incident Dementia in GenScot and LBC1936") + theme(plot.title = element_text(hjust = 0.5, size = 10))+
  #guides(color = "none", size = "none")+
  #guides(shape = guide_legend(override.aes = list(size=5))) + 
  coord_flip()
dev.off()


