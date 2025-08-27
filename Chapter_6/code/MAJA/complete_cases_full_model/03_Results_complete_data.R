#load dyplr and data.table packages 
library(dplyr)
library(data.table)

#read in betas 
mean_beta <- read.csv("mean_beta.csv")

#check data dimensions 
dim(mean_beta)
#[1] 439   6

#read in protein data to use column names as ids 
prot_id <- readRDS("GS_MS_prots_MAJA_cog_adPRS_cogPRS_fully_corrected.rds")

#reintroduce the column names for outcomes 
names(mean_beta) <- c( "digit symbol", "verbal fluency", "vocabulary","logical memory", "AD PRS", "cognitive PRS")

#reintroduce the protein names as rownames 
mean_beta$Protein <- colnames(prot_id) 

#read in PIPs (posterior inclusion probabilities)
pip <- read.table("mean_prob.txt")

#check data dimensions 
dim(pip)
#[1] 439   1

#name column 'PIP'
names(pip)[1] <- "PIP"

#bind together betas and pips
data <- cbind(mean_beta,pip)


#load tidyr package
library(tidyr)

#pivot dataset longer to merge in variance for each beta 
data <- data %>% 
  pivot_longer(
    cols = c(1:6),
    names_to = "Trait", 
    values_to = "Beta")

#read in the variance of the betas 
var <- read.csv("var_beta.csv")

#check data dimensions 
dim(var)
#[1] 439   6

#set columns names to outcomes 
names(var) <- c( "digit_symbol", "verbal total", "vocabulary","LM", "AD_PRS", "cog_PRS")

#get the square root of the variances to calculate standard deviation 
var$digit_symbol <- sqrt(var$digit_symbol)
var$verbal_total <- sqrt(var$verbal_total)
var$vocabulary <- sqrt(var$vocabulary)
var$LM <- sqrt(var$LM)
var$AD_PRS <- sqrt(var$AD_PRS)
var$cog_PRS <- sqrt(var$cog_PRS)

#make a protein id column 
var$Protein <- colnames(prot_id)

#rename column names for plotting 
names(var) <- c( "digit symbol", "verbal fluency", "vocabulary","logical memory", "AD PRS", "cognitive PRS", "Protein")

#pivot dataset to merge with betas/pips 
var <- var %>% 
  pivot_longer(
    cols = c(1:6),
    names_to = "Trait", 
    values_to = "sd")


#join SDs dataset with betas/pips by trait and protein 
joined_df <- left_join(data, var, by=c('Trait', 'Protein'))

#calculate upper and lower SD limits
joined_df$UpperSD <- joined_df$Beta + joined_df$sd 
joined_df$LowerSD <- joined_df$Beta - joined_df$sd 


#read in protein annotation file
ids <- read.csv("short_annots_final.csv")

#set id column name to protein 
names(ids)[1] <- "Protein"

#merge in the annotations 
joined_df <- left_join(joined_df, ids, by = "Protein")



joined_df$Significance <- ifelse(joined_df$PIP >= 0.95 & (joined_df$Beta > 0 & joined_df$LowerSD > 0 | joined_df$Beta < 0 & joined_df$UpperSD < 0), "Significant", "Non significant")

table(joined_df$Significance)

# Non significant     Significant
# 2607              27


#filter to rows with PIP greater than or equal to 95%
sig <- joined_df[joined_df$PIP >= 0.95, ]

#check number of sig hits based on pip 
dim(sig)
#[1] 42  9

table(sig$Significance)

# Non significant     Significant
# 15              27


#check number of unqiue proteins 
length(unique((sig$Name)))
#[1] 7


#save out results 
write.csv(joined_df, "GS_four_cog_tests_ADPRS_cogPRS_MAJA_all_complete_data_results.csv") 
write.csv(sig, "/GS_four_cog_tests_ADPRS_cogPRS_MAJA_sig_complete_data_results.csv") 




library(ggplot2)
library(dplyr)
library(tidytext)  # For reorder_within and scale_y_reordered

# Prepare data
sig <- sig %>%
  mutate(
    trait_label = Trait,  # original Trait used for color
    Trait = reorder_within(Trait, Beta, Name)  # reordered version used for y-axis
  )

# Define Okabe-Ito palette (6 outcomes)
okabe_ito_palette_no_yellow <- c(
  "#E69F00", "#56B4E9", "#009E73",
  "#0072B2", "#D55E00", "#CC79A7"
)

# Create TIFF
setwd("")
tiff("multi_panel_hits.tiff", width = 14 * 300, height = 12 * 300, res = 300)

ggplot(sig, aes(x = Beta, y = Trait, color = trait_label, shape = Significance)) +
  geom_point(size = 3) +
  geom_errorbar(aes(xmax = UpperSD, xmin = LowerSD), size = 1, width = 0.2) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
  facet_wrap(~ Name, scales = "free_y", ncol = 3, labeller = label_wrap_gen(width = 20)) +  # Wrap long facet titles
  scale_y_reordered() +
  scale_color_manual(values = okabe_ito_palette_no_yellow, guide = "none") +  # Remove trait legend
  scale_shape_manual(
    values = c("Significant" = 17, "Non significant" = 16),
    name = "Significance"
  ) +
  labs(x = "Beta [+/- SD]", y = "Outcome") +
  theme_minimal(base_size = 16) +
  theme(
    legend.position = "right",
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_line(color = "gray70", size = 0.5),
    plot.margin = margin(t = 10, r = 30, b = 10, l = 10)
  )

dev.off()
