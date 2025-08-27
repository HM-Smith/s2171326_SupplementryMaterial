#read in results 
full = read.csv("Dementia_GS_MS_prots_KNNimpute_model_nominally_significant.csv")

complete = read.csv("Dementia_GS_MS_prots_full_model_nominally_significant.csv")

#make uniquw list of proteins from both datasets
sig_proteins = union(complete$Name, full$Name)

length(sig_proteins)
#[1] 47

#read in full results datasets 
full_all = read.csv("Dementia_GS_MS_prots_KNNimpute_model_supplementary.csv")
complete_all = read.csv("Dementia_GS_MS_prots_full_model_supplementary.csv")

#subset to unique list
df_complete_sig <- complete_all[complete_all$Name %in% sig_proteins, ]
df_full_sig <- full_all[full_all$Name %in% sig_proteins, ]

dim(df_complete_sig)
#[1] 47 12
dim(df_full_sig)
#[1] 47 12

#make model identifier
df_complete_sig$Model = rep("sensitivity", 47)
df_full_sig$Model = rep("full", 47)

#bind datasets together
df = rbind(df_complete_sig, df_full_sig)

#set significance threshold
bonf = 0.05/161 

#make significance identfier 
df$Significance = ifelse(df$P < bonf, "bonferroni significant", ifelse(df$P < 0.05, "nominally significant", "non significant"))

#order by effect size
df$Name = factor(df$Name, levels=unique(df$Name[order(df$Hazard_Ratio)]))

library(ggplot2)
#plot time-to-dementia results 
setwd("")
tiff("Dementia_complete_full_compar.tiff", width = 8*300, height = 8*300, res = 300)

pd <- position_dodge(width = 0.6)

ggplot(df, aes(y = Hazard_Ratio, x = Name, color = Model, shape = Significance)) +
  geom_point(size = 2.5, position = pd) +
  geom_errorbar(aes(ymin = ci.lower, ymax = ci.upper), size = 0.5, width = 0, position = pd) +
  scale_color_manual(values = c("#E69F00", "#56B4E9")) +
  scale_shape_manual(values = c(17, 16, 0)) +
  geom_hline(yintercept = 1, size = 0.25, linetype = "dashed") +
  coord_flip() +
  theme_classic() +
  ylab("Hazard Ratio [95% CI]") +
  xlab("Protein") +
  guides(shape = guide_legend(override.aes = list(size = 5)))

dev.off()