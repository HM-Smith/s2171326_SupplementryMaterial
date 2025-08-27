#read in results
basic = read.csv("Dementia_GS_MS_prots_basic_model_nominally_significant.csv")

full = read.csv("/Dementia_GS_MS_prots_KNNimpute_model_nominally_significant.csv")

#unique list of proteins in from both datasets
sig_proteins = union(basic$Name, full$Name)

#check no of proteins
length(sig_proteins)
#[1] 34

basic_all = read.csv("Dementia_GS_MS_prots_basic_model_supplementary.csv")

full_all = read.csv("Dementia_GS_MS_prots_KNNimpute_model_supplementary.csv")

#subset to proteins in unique list
df_basic_sig <- basic_all[basic_all$Name %in% sig_proteins, ]
df_full_sig <- full_all[full_all$Name %in% sig_proteins, ]

dim(df_basic_sig)
#[1] 34 12
dim(df_full_sig)
#[1] 34 12

#add model identifier
df_basic_sig$Model = rep("basic", 34)
df_full_sig$Model = rep("full", 34)

#bind datasets together 
df = rbind(df_basic_sig, df_full_sig)

#set sig threshold 
bonf = 0.05/161 

#make significance identifier 
df$Significance = ifelse(df$P < bonf, "bonferroni significant", ifelse(df$P < 0.05, "nominally significant", "non significant"))

#order by effect size
df$Name = factor(df$Name, levels=unique(df$Name[order(df$Hazard_Ratio)]))

library(ggplot2)
#plot time-to-dementia results 
setwd("")
tiff("Dementia_basic_full_compar.tiff", width = 8*300, height = 8*300, res = 300)

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