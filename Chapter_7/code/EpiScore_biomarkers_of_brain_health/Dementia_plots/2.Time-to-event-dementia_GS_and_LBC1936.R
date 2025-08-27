library(tidyverse)

#load FDR significant results for GS
GenScot_FDR <- read.csv("Coxme_Dementia_GenScot_84_episcores_FDR_sig_basic_16082023.csv")

#Check dimensions 
dim(GenScot_FDR)
#[1] 13 18

#load nominally significant results for LBC1936
LBC1936_Nom_sig <- read.csv("coxPH_Dementia_LBC1936_84_EpiScores_basic_model_nominally_sig_20072023.csv")

#subset LBC1936 nominal results to FDR significant scores in GS
LBC1936_Nom_sig_common <- LBC1936_Nom_sig[LBC1936_Nom_sig$SeqId %in% GenScot_FDR$SeqId, ]

#check how many scores present
dim(LBC1936_Nom_sig_common)
#[1]  0 18

#load all LBC1936 results 
LBC1936_all <- read.csv("coxPH_Dementia_LBC1936_84_EpiScores_basic_model_WITH_IDS_20072023.csv")

#subset to FDR significant results in GS 
LBC1936_common <- LBC1936_all[LBC1936_all$SeqId %in% GenScot_FDR$SeqId,]

#check that all are present 
dim(LBC1936_common)
#[1]  13 18

#add cohort and significance columns 
LBC1936_common$Cohort <- rep("LBC1936", 13)
LBC1936_common$Significance <- rep("Non significant", 13)

#add cohort and significance columns 
GenScot_FDR$Cohort <- rep("GS", 13)
GenScot_FDR$Significance <- rep("FDR significant", 13)


#select columns of interest
GenScot_FDR <- GenScot_FDR %>% select(SeqId, Hazard_Ratio, ci.lower, ci.upper, Gene_Name, Cohort, Significance, Model)
LBC1936_common <- LBC1936_common %>% select(SeqId, Hazard_Ratio, ci.lower, ci.upper, Gene_Name, Cohort, Significance, Model)

#bind results together
all<- rbind(GenScot_FDR, LBC1936_common)

#order by effect size 
all$Gene_Name = factor(all$Gene_Name, levels=unique(all$Gene_Name[order(all$Hazard_Ratio)]))

#plot time-to-dementia results 
setwd("")
pdf("time_to_event_All_3_cohorts_in_relation_to_FDR_sig_GenScot_16082023.pdf")

  ggplot(all, aes(y = Hazard_Ratio, x = Gene_Name, color = Cohort, shape = Significance)) +
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
  #ggtitle("Time-to-dementia in GS and LBC1936") + theme(plot.title = element_text(hjust = 0.5, size = 10))+
  #guides(color = "none", size = "none")+
  guides(shape = guide_legend(override.aes = list(size=5))) + 
  coord_flip()

dev.off()




