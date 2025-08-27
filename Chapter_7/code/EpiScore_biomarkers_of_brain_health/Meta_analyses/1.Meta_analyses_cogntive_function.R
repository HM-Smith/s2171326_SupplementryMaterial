library(metafor)
library(dplyr)

#read in cognitive results for each cohort
genscot_cog <- read.csv("GS20K_G_intercept_assocs_84_EpiScores_with_ids_basic_20012023.csv")

LBC1921_cog <- read.csv("G_intercept_No_Domain_LBC1921_assocs_84_EpiScores_with_ids_basic_model_21042023.csv")

LBC1936_cog <- read.csv("G_intercept_No_Domain_LBC1936_assocs_84_EpiScores_with_ids_basic_20072023.csv")

#select columns of interest
genscot_cog <- genscot_cog %>% select(4:20)
LBC1921_cog <- LBC1921_cog %>% select(3:19)
LBC1936_cog <- LBC1936_cog %>% select(3:19)

#rename columns to match
names(genscot_cog)[13] <- "Gene_Name"
names(LBC1936_cog)[13] <- "Gene_Name"

#make list of episcores
EpiScores <- genscot_cog$SeqId 

#set rownames as episcores for each results file
rownames(genscot_cog) <- genscot_cog$SeqId
rownames(LBC1921_cog) <- LBC1921_cog$SeqId
rownames(LBC1936_cog) <- LBC1936_cog$SeqId

#create empty results file for meta-analyses
outcomes <- data.frame("SeqId" = EpiScores, "Gene_Name" = NA, "Target" = NA, "Uniprot"= NA, "q" = NA, "p_for_q" = NA, "I^2" = NA, "tau^2" = NA, "setau^2" = NA, "beta" = NA, 
                       "se" = NA, "z" = NA, "p" = NA, "H^2" = NA, "lower_CI" = NA, "Upper_CI" = NA)

#set episcore names as rownames
rownames(outcomes) <- outcomes$SeqId


# for loop for looping over episcores and performing meta-analysis 
for (i in EpiScores) {
  loopdata <- rbind(genscot_cog[i,], LBC1921_cog[i,], LBC1936_cog[i,]) # a dataframe with the episcores per row and columns with betas and SE, and can loop over the rows
  
  
  res <- rma(beta, sei = SE, data = loopdata) # meta-analysis
  
  outcomes[i,1] <- loopdata$SeqId[1]
  outcomes[i,2] <- loopdata$Gene_Name[1]
  outcomes[i,3] <- loopdata$Target[1]
  outcomes[i,4] <- loopdata$UniProt[1]
  
  outcomes[i,5] <- as.numeric(res[17]) # q statistic 
  outcomes[i,6] <- as.numeric(res[18]) # p value for q
  outcomes[i,7] <- as.numeric(res[13]) # I^2
  outcomes[i,8] <- as.numeric(res[9]) # tau^2
  outcomes[i,9] <- as.numeric(res[10]) # setau^2
  outcomes[i,10] <- as.numeric(res[2]) #beta
  outcomes[i,11] <- as.numeric(res[3]) # se
  outcomes[i,12] <- as.numeric(res[4]) # z
  outcomes[i,13] <- as.numeric(res[5]) # p value
  outcomes[i,14] <- as.numeric(res[14]) # H^2
  outcomes[i,15] <- as.numeric(res[6])# lower CI
  outcomes[i,16] <- as.numeric(res[7])# upper CI 
  
  }

#create FDR corrected column
outcomes$FDR <- p.adjust(outcomes$p, method = "BH")

#save out results
write.csv(outcomes, file = "cognitive_ability_all_three_cohorts_meta_analysis_basic_results_20072023.csv")

#filter to FDR significant results
outcomes_sig <- outcomes %>% filter(FDR < 0.05)

#check number of FDR significant results
dim(outcomes_sig)
#[1] 36 17

#save out significant results
write.csv(outcomes_sig, file = "cognitive_ability_all_three_cohorts_meta_analysis_significant_results_basic_models_20072023.csv")

#order by effect size
outcomes_sig$Gene_Name = factor(outcomes_sig$Gene_Name, levels=unique(outcomes_sig$Gene_Name[order(outcomes_sig$beta)]))


#plot meta-analysis results 
setwd("")
pdf("Cognitive_ability(intercept)_meta_analysis_all_3_cohorts_basic_models_ordered_beta_21072023.pdf")
library(ggplot2)

ggplot(outcomes_sig, aes(y = beta, x = Gene_Name)) +
  geom_point(size = 2.5, position=position_dodge(0.5), colour = "#F8766D") +
  geom_errorbar(aes(ymax = Upper_CI, ymin = lower_CI), size = .6, width =
                  0, position=position_dodge(0.5), colour = "#F8766D") +
  geom_hline(aes(yintercept = 0), size = .25, linetype = "dashed") +
  #coord_trans(y = scales:::exp_trans(10)) + 
  theme_classic() +
  #theme(panel.grid.minor = element_blank()) +
  ylab("Beta [95% CI]") +
  xlab("EpiScore") +
  ggtitle("Meta-analysis of Cognitive ability in GS, LBC1921,and LBC1936") +
  theme(plot.title = element_text(hjust = 0.5, size = 12))+
  #guides(color = "none", size = "none")+
  #guides(shape = guide_legend(override.aes = list(size=5))) + 
  coord_flip()
dev.off()
###################################################################################################################
###################################################################################################################

library(metafor)
library(dplyr)

#read in cognitive results for all cohorts - full models 
genscot_cog <- read.csv("GS20K_G_intercept_assocs_84_EpiScores_with_ids_full_model_no_yrsedu_31032023.csv")

LBC1921_cog <- read.csv("G_intercept_No_Domain_LBC1921_assocs_84_EpiScores_with_ids_full_model_no_yrsedu_21042023.csv")

LBC1936_cog <- read.csv("G_intercept_No_Domain_LBC1936_assocs_84_EpiScores_with_ids_full_model_no_yrsedu_27102023.csv")

#select columns of interest
genscot_cog <- genscot_cog %>% select(4:20)
LBC1921_cog <- LBC1921_cog %>% select(3:19)
LBC1936_cog <- LBC1936_cog %>% select(3:19)

#rename columns to match 
names(genscot_cog)[13] <- "Gene_Name"
names(LBC1936_cog)[13] <- "Gene_Name"

#make episcore list
EpiScores <- genscot_cog$SeqId 

#set rownames for each of the results files to episcores
rownames(genscot_cog) <- genscot_cog$SeqId
rownames(LBC1921_cog) <- LBC1921_cog$SeqId
rownames(LBC1936_cog) <- LBC1936_cog$SeqId
z

#make empty dataframe for meta-analysis results 
outcomes <- data.frame("SeqId" = EpiScores, "Gene_Name" = NA, "Target" = NA, "Uniprot"= NA, "q" = NA, "p_for_q" = NA, "I^2" = NA, "tau^2" = NA, "setau^2" = NA, "beta" = NA, 
                       "se" = NA, "z" = NA, "p" = NA, "H^2" = NA, "lower_CI" = NA, "Upper_CI" = NA)

#set episcores as rownames
rownames(outcomes) <- outcomes$SeqId

# for loop for looping over episcores and performing meta-analysis 
for (i in EpiScores) {
  loopdata <- rbind(genscot_cog[i,], LBC1921_cog[i,], LBC1936_cog[i,]) #dataframe with the episcores per row and columns with betas and variances, and can loop over the rows
  
  
  res <- rma(beta, sei = SE, data = loopdata) #meta-analysis
  
  outcomes[i,1] <- loopdata$SeqId[1]
  outcomes[i,2] <- loopdata$Gene_Name[1]
  outcomes[i,3] <- loopdata$Target[1]
  outcomes[i,4] <- loopdata$UniProt[1]
  
  outcomes[i,5] <- as.numeric(res[17]) # q statistic
  outcomes[i,6] <- as.numeric(res[18]) # p value for q
  outcomes[i,7] <- as.numeric(res[13]) # I^2
  outcomes[i,8] <- as.numeric(res[9]) # tau^2
  outcomes[i,9] <- as.numeric(res[10]) # setau^2
  outcomes[i,10] <- as.numeric(res[2]) #beta
  outcomes[i,11] <- as.numeric(res[3]) # se
  outcomes[i,12] <- as.numeric(res[4]) # z
  outcomes[i,13] <- as.numeric(res[5]) # p value
  outcomes[i,14] <- as.numeric(res[14]) # H^2
  outcomes[i,15] <- as.numeric(res[6])# lower CI
  outcomes[i,16] <- as.numeric(res[7])# upper CI 
  
}

#create FDR corrected column 
outcomes$FDR <- p.adjust(outcomes$p, method = "BH")

#save out all results 
write.csv(outcomes, file = "cognitive_ability_all_three_cohorts_meta_analysis_full_no_years_edu_results_30102023.csv")

#filter to FDR significant results 
outcomes_sig <- outcomes %>% filter(FDR < 0.05)

#check number of FDR significant results
dim(outcomes_sig)
#[1] 18 17

#save out FDR significant results
write.csv(outcomes_sig, file = "cognitive_ability_all_three_cohorts_meta_analysis_significant_results_full_no_years_edu_models_30102023.csv")

#order by effect size
outcomes_sig$Gene_Name = factor(outcomes_sig$Gene_Name, levels=unique(outcomes_sig$Gene_Name[order(outcomes_sig$beta)]))

#plot meta-analyses - full model
setwd("")
pdf("Cognitive_ability(intercept)_meta_analysis_all_3_cohorts_full_no_years_edu_models_ordered_betas_meta_only_30102023.pdf")
library(ggplot2)

ggplot(outcomes_sig, aes(y = beta, x = Gene_Name)) +
  geom_point(size = 2.5, position=position_dodge(0.5), colour = "#F8766D") +
  geom_errorbar(aes(ymax = Upper_CI, ymin = lower_CI), size = .6, width =
                  0, position=position_dodge(0.5), colour = "#F8766D") +
  geom_hline(aes(yintercept = 0), size = .25, linetype = "dashed") +
  theme_classic() +
  ylab("Beta [95% CI]") +
  xlab("EpiScore") +
  ggtitle("Meta-analysis of Cognitive ability in GS, LBC1921,and LBC1936") +
  theme(plot.title = element_text(hjust = 0.5, size = 12))+
  expand_limits(y = 0.1) +
  coord_flip()
dev.off()
