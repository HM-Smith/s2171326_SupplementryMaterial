# Load necessary library
library(dplyr)

#load phenotype file 
load("20191125_1614_pheno_combined.Rdata")

#remove studies that are not required
pheno_combined = pheno_combined[!(pheno_combined$STUDY %in% c("Harmony","IPT3", "Parkinson")), ]

#Load outcome
response = read.csv("")

#subset to complete data
response = response[complete.cases(response),]

#remove studies that are not required
response = response[!(response$study %in% c("SATSA7")), ]

#set id as character
response$twinnr = as.character(response$twinnr)

#rename studies 
response$study <- recode(response$study,
                         "SATSA10" = "IPT10",
                         "SATSA5"  = "IPT5",
                         "SATSA6"  = "IPT6",
                         "SATSA8"  = "IPT8",
                         "SATSA9"  = "IPT9")

#rename columns 
names(response)[2] <- "TWINNR"
names(response)[6] <- "STUDY"

#join datasets
response = left_join(response, pheno_combined, by = c("TWINNR", "STUDY"))

#Load methylation data
load("20180128_0046_betas_ct_batch.Rdata")

#transpose meth 
betas_ct_batch = t(betas_ct_batch)

#subset protein data to those with meth
response = response[response$Sample %in% rownames(betas_ct_batch),]

dim(response)
#[1] 896  18

#Get unique twin pair IDs
unique_pairs <- response %>%
  select(pairid) %>%
  distinct()

dim(unique_pairs)



# Randomly sample 80% of the pairs for training
# Split the original data
set.seed(123)  # for reproducibility
train_pair_ids <- sample(unique_pairs$pairid, size = 0.8 * nrow(unique_pairs))

train_data <- response %>% filter(pairid %in% train_pair_ids)
test_data <- response %>% filter(!pairid %in% train_pair_ids)

# Check results
cat("Training samples:", nrow(train_data), "\n")
cat("Testing samples:", nrow(test_data), "\n")

# Check results
cat("Training samples ids:", length(unique(train_data$TWINNR)), "\n")
cat("Testing samples ids:",length(unique(test_data$TWINNR)), "\n")


#descriptives of training and test sets 

summary(train_data$bio_age)

summary(train_data$AGE)


summary(test_data$bio_age)

  
summary(test_data$AGE)

table(train_data$SEX)


table(test_data$SEX)



table(train_data$SMOKE)


table(test_data$SMOKE)



write.csv(train_data, file = "")
write.csv(test_data, file = "")

