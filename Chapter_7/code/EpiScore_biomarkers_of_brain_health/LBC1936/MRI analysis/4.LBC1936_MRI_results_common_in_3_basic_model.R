library(readr)
library(tidyverse)

#read in results files
Brain <- read.csv("Totalbrain_intercept_LBC1936_assocs_84_EpiScores_basic_15092023.csv")
GM <- read.csv("GreyMatter_intercept_LBC1936_assocs_84_EpiScores_basic_15092023.csv")
NAWM <- read.csv("NAWM_intercept_LBC1936_assocs_84_EpiScores_basic_15092023.csv")
WMH <- read.csv("WMH_intercept_LBC1936_assocs_84_EpiScores_basic_15092023.csv")


#edit episcore names for merging with episcore ids 
Brain$SeqId <- gsub("X", "", as.character(Brain$SeqId))
Brain$SeqId <- gsub("\\.", "-", as.character(Brain$SeqId))

GM$SeqId <- gsub("X", "", as.character(GM$SeqId))
GM$SeqId <- gsub("\\.", "-", as.character(GM$SeqId))

NAWM$SeqId <- gsub("X", "", as.character(NAWM$SeqId))
NAWM$SeqId <- gsub("\\.", "-", as.character(NAWM$SeqId))

WMH$SeqId <- gsub("X", "", as.character(WMH$SeqId))
WMH$SeqId <- gsub("\\.", "-", as.character(WMH$SeqId))

#read in episcore ids 
library(readxl)
Epi_ID <-read_excel("Index_109.xlsx")
Epi_ID <- Epi_ID %>% rename(SeqId = `Identifier (SOMAScan SeqId or Olink name)`,
                            Gene_Name = 'Gene Name')

#merge with episcore ids
Brain <- Brain %>% left_join(Epi_ID, by = "SeqId")
GM <- GM %>% left_join(Epi_ID, by = "SeqId")
NAWM <- NAWM %>% left_join(Epi_ID, by = "SeqId")
WMH<- WMH %>% left_join(Epi_ID, by = "SeqId")

#make FDR corrected column
Brain$FDR <- p.adjust(Brain$P, method = 'BH')
GM$FDR <- p.adjust(GM$P, method = 'BH')
NAWM$FDR <- p.adjust(NAWM$P, method = 'BH')
WMH$FDR <- p.adjust(WMH$P, method = 'BH')

#make measure column
Brain$Measure <- rep("Total Brain Volume", 84)
GM$Measure <- rep("Grey matter Volume", 84)
NAWM$Measure <- rep("Normal appearing White matter Volume", 84)
WMH$Measure <- rep("White matter hyperintensity volume", 84)

#filter to FDR significant results
Brain_sig <- Brain %>% filter(FDR < 0.05)
GM_sig <- GM %>% filter(FDR < 0.05)
NAWM_sig <- NAWM %>% filter(FDR < 0.05)
WMH_sig <- WMH %>% filter(FDR < 0.05)


#check number of FDR significant results
dim(Brain_sig)
#[1] 21 20
dim(GM_sig)
#[1] 28 20
dim(NAWM_sig)
#[1] 17 20
dim(WMH_sig)
#[1]  3 20


#filter to nominally significant results 
Brain_nom <- Brain %>% filter(P < 0.05)
GM_nom <- GM %>% filter(P < 0.05)
NAWM_nom <- NAWM %>% filter(P < 0.05)
WMH_nom <- WMH %>% filter(P < 0.05)

#check number of nominally significant results 
dim(Brain_nom)
#[1] 35 20
dim(GM_nom)
#[1] 35 20
dim(NAWM_nom)
#[1] 24 20
dim(WMH_nom)
#[1]  9 20


#check how many episcores associate with 3 or more measures
check <- plyr::count(c(Brain_sig$SeqId,GM_sig$SeqId, NAWM_sig$SeqId, WMH_sig$SeqId))
check <- check %>% filter(freq > 2) #11

#filter to episcores that associate with 3 or more measures
Brain_sig_c <- Brain_sig[Brain_sig$SeqId %in% check$x, ]#10
GM_sig_c <- GM_sig[GM_sig$SeqId %in% check$x, ]#11
NAWM_sig_c <- NAWM_sig[NAWM_sig$SeqId %in% check$x, ]#11
WMH_sig_c <- WMH_sig[WMH_sig$SeqId %in% check$x, ]#3

#bind the filtered files together
data <- rbind(Brain_sig_c, GM_sig_c)
data <- rbind(data, NAWM_sig_c)
data <- rbind(data, WMH_sig_c)


#postively recode WMH volume 
data$beta_new <- ifelse(data$phenotype == "IWMHV" & data$beta > 0, -data$beta, ifelse(data$phenotype == "IWMHV" & data$beta < 0, abs(data$beta), data$beta))
data$CI_upper_new <- ifelse(data$phenotype == "IWMHV" & data$beta > 0, -data$ci.upper, ifelse(data$phenotype == "IWMHV" & data$beta < 0, abs(data$ci.upper), data$ci.upper))
data$CI_lower_new <- ifelse(data$phenotype == "IWMHV" & data$beta > 0, -data$ci.lower, ifelse(data$phenotype == "IWMHV" & data$beta < 0, abs(data$ci.lower), data$ci.lower))


#order by effect size
data$Gene_Name = factor(data$Gene_Name, levels=unique(data$Gene_Name[order(data$beta_new)]))


#plot episcores that associate with 3 or more measures
setwd("")
pdf("Baseline_MRI_basic_models_common_atleast_3_beta_recoded_20092023.pdf")
ggplot(data, aes(y = beta_new, x = Gene_Name, color = Measure)) +
  geom_point(size = 2.5, position=position_dodge(0.7)) +
  geom_errorbar(aes(ymax = CI_upper_new, ymin = CI_lower_new), size = .6, width =
                  0, position=position_dodge(0.7)) +
  geom_hline(aes(yintercept = 0), size = .25, linetype = "dashed") +
  #coord_trans(y = scales:::exp_trans(10)) + 
  theme_classic() +
  #theme(panel.grid.minor = element_blank()) +
  ylab("Beta [95% CI]") +
  xlab("EpiScore") +
  ggtitle("Baseline MRI measures of brain health in LBC1936") +
  theme(plot.title = element_text(hjust = 0.5, size = 12))+
  #guides(color = "none", size = "none")+
  #guides(shape = guide_legend(override.aes = list(size=5))) + 
  coord_flip()
dev.off()


