library(tidyverse)
library(ggplot2)

#load FDR significant results for LBC1921
LBC1921_FDR <- read.csv("LBC1921_dementia_84_episcores_FDR_sig_logistic_basic_model_05072023.csv")

# Add cohort and significance columns 
LBC1921_FDR$Cohort <- rep("LBC1921", 3)
LBC1921_FDR$Significance <- rep("FDR significant", 3)

#load LBC1936 nominally significant results
LBC1936_Nom_sig <- read.csv("logistic_Dementia_LBC1936_84_EpiScores_basic_model_nominally_sig_20072023.csv")

#subset to EpiScores in LBC1921 FDR 
LBC1936_Nom_sig_common <- LBC1936_Nom_sig[LBC1936_Nom_sig$SeqId %in% LBC1921_FDR$SeqId, ]

#check if any EpiScores nominally significant in LBC1936 that were FDR significant in LBC1921
dim(LBC1936_Nom_sig_common)
#[1]  0 16

#load all results for LBC1936
LBC1936_all <- read.csv("logistic_Dementia_LBC1936_84_EpiScores_basic_model_WITH_IDS_20072023.csv")

#subset to FDR significant scores in LBC1921
LBC1936_common <- LBC1936_all[LBC1936_all$SeqId %in% LBC1921_FDR$SeqId,]

#check all scores are present
dim(LBC1936_common)
#[1]  3 16

#add cohort and significance columns 
LBC1936_common$Cohort <- rep("LBC1936", 3)
LBC1936_common$Significance <- rep("Non significant", 3)


#load GS nominally significant results 
GenScot_Nom_sig<- read.csv("logistic_Dementia_GenScot_84_episcores_nominally_sig_basic_16082023.csv")

#subset to FDR significant results in LBC1921 
GenScot_Nom_sig <- GenScot_Nom_sig[GenScot_Nom_sig$SeqId %in% LBC1921_FDR$SeqId, ]

#check if any scores are nominally significant
dim(GenScot_Nom_sig)
#[1]  1 16

#add cohort and significance columns 
GenScot_Nom_sig$Cohort <- "GS"
GenScot_Nom_sig$Significance <- "Nominally significant"

#load to all GS results
GenScot_all <- read.csv("logistic_Dementia_GenScot_84_episcores_with_ids_basic_16082023.csv")

#subset to FDR significant results in LBC1921 but not in GS nominally significant results
GenScot_all_common <- GenScot_all[GenScot_all$SeqId %in% LBC1921_FDR$SeqId & !GenScot_all$SeqId %in% GenScot_Nom_sig$SeqId,]

#check scores are present 
dim(GenScot_all_common)
#[1]  2 16

#add cohort and significance columns 
GenScot_all_common$Cohort <- rep("GS", 2)
GenScot_all_common$Significance <- rep("Non significant", 2)

#rename to Odds Ratio in GS
names(GenScot_all_common)[4] <- "Odds_Ratio"
names(GenScot_Nom_sig)[4] <- "Odds_Ratio"

#select columns of interest 
LBC1921_FDR <- LBC1921_FDR %>% select(SeqId, Odds_Ratio, ci.lower, ci.upper, Gene_Name, Target, Cohort, Significance)
LBC1936_common <- LBC1936_common %>% select(SeqId, Odds_Ratio, ci.lower, ci.upper, Gene_Name, Target, Cohort, Significance)
GenScot_all_common <- GenScot_all_common %>% select(SeqId, Odds_Ratio, ci.lower, ci.upper, Gene_Name, Target, Cohort, Significance)
GenScot_Nom_sig <- GenScot_Nom_sig %>% select(SeqId, Odds_Ratio, ci.lower, ci.upper, Gene_Name, Target, Cohort, Significance)

#bind results together
all<- rbind(LBC1921_FDR, LBC1936_common)
all <- rbind(all, GenScot_all_common)
all <- rbind(all, GenScot_Nom_sig)

#order by effect size 
all$Gene_Name = factor(all$Gene_Name, levels=unique(all$Gene_Name[order(all$Odds_Ratio)]))


# plot logistic regression results with three cohorts on the same plot
setwd("")
pdf("Logistic_All_3_cohorts_in_relation_to_FDR_sig_LBC1921_16082023.pdf")
ggplot(all, aes(y = Odds_Ratio, x = Gene_Name, color = Cohort, shape = Significance)) +
  geom_point(size = 2.5, position=position_dodge(0.5)) +
  geom_errorbar(aes(ymax = ci.upper, ymin = ci.lower),size = .5, width =
                  0, position=position_dodge(0.5)) +
  scale_shape_manual(values = c(17,16,0)) +
  geom_hline(aes(yintercept = 1), size = .25, linetype = "dashed") +
  #coord_trans(y = scales:::exp_trans(10)) + 
  theme_classic() +
  #theme(panel.grid.minor = element_blank()) +
  ylab("Odds Ratio [95% CI]") +
  xlab("EpiScore") +
  #ggtitle("Dementia in GS, LBC1936 and LBC1921") + theme(plot.title = element_text(hjust = 0.5, size = 10))+
  #guides(color = "none", size = "none")+
  guides(shape = guide_legend(override.aes = list(size=5))) + 
  coord_flip()
dev.off()


