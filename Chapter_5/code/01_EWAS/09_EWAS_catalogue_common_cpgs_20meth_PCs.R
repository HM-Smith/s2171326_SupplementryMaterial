library(dplyr)

#read in common cpgs
df_meth <- read.csv("common_CpGs_20methPCs_18411_full_models.csv")

#check dimensions 
dim(df_meth)
#[1] 36  4

#make list of cpg names
cpgs<- df_meth$Probe

#check its same length
length(cpgs)
#[1] 36

#read in BMI EWAS cat studies
bmi_cat <- read.csv("BMI_EWAS_cat.csv")

#check data dimensions
dim(bmi_cat)
#[1] 6427   33

#filter cat to whole blood, sample size > 1000 and P < 3.8e-8
bmi_cat <- bmi_cat[bmi_cat$Tissue == "Whole blood" & bmi_cat$P < 3.6e-8 & bmi_cat$N > 1000, ]

#check data dimensions 
dim(bmi_cat)
#[1] 1388   33

#check unique study ids 
unique(bmi_cat$PMID)
#[1] 25935004 28002404 28095459 26110892 26119815 31510868 29278407 29099282
#[9] 33550919 33239103

#check for overlaping assocs between our study and EWAS cat
bmi_in_cat <- bmi_cat[bmi_cat$CpG %in% cpgs,]

#check data dimensions
dim(bmi_in_cat)
#[1] 86 33

#check num of unique cpgs that overlap
length(unique(bmi_in_cat$CpG))
#[1] 11

#print the cpg names
unique(bmi_in_cat$CpG)
#[1] "cg14476101" "cg16246545" "cg17901584" "cg23032421" "cg10975897"
#[6] "cg00574958" "cg17058475" "cg11024682" "cg06500161" "cg08309687"
#[11] "cg27243685"


#print the unique pubmed ids 
unique(bmi_in_cat$PMID)
#[1] 28095459 28002404 26110892 25935004 26119815 31510868 29278407 29099282



#read in EWAS cat with glucose studies
Glucose_cat <- read.csv("Glucose_ewas_cat.csv")

#check data dimensions 
dim(Glucose_cat)
#[1] 3162   33

#filter cat to whole blood, sample size > 1000 and P < 3.8e-8
Glucose_cat <- Glucose_cat[Glucose_cat$Tissue == "Whole blood" & Glucose_cat$P < 3.6e-8 & Glucose_cat$N > 1000, ]

#check data dimensions 
dim(Glucose_cat)
#[1]  3 33

#check unique pubmed ids
unique(Glucose_cat$PMID)
#[1] 27019061

#check overlap between our study and EWAS cat 
Glucose_in_cat <- Glucose_cat[Glucose_cat$CpG %in% cpgs,]

#check data dimensions 
dim(Glucose_in_cat)
#[1]  3 33

#check unique cpgs that overlap
length(unique(Glucose_in_cat$CpG))
#[1] 2

#print unique cpg names
unique(Glucose_in_cat$CpG)
#[1] "cg11024682" "cg06500161"

#print unique pubmed id
unique(Glucose_in_cat$PMID)
#[1] 27019061

#read in EWAS cat with HDL cholesterol studies
HDL_cat <- read.csv("HDL_ewas_cat.csv")

#check data dimensions 
dim(HDL_cat)
#[1] 1026   33

#filter cat to whole blood, sample size > 1000 and P < 3.8e-8
HDL_cat <- HDL_cat[HDL_cat$Tissue == "Whole blood" & HDL_cat$P < 3.6e-8 & HDL_cat$N > 1000, ]

#check data dimensions
dim(HDL_cat)
#[1] 43 33

#print unique pubmed ids
unique(HDL_cat$PMID)
#[1] 28194238 25583993

#check overlap between our study and EWAS cat
HDL_in_cat <- HDL_cat[HDL_cat$CpG %in% cpgs,]

#check data dimensions 
dim(HDL_in_cat)
#[1]  8 33

#check num of unique cpgs
length(unique(HDL_in_cat$CpG))
#[1] 4

#print unique cpg names
unique(HDL_in_cat$CpG)
#[1] "cg17901584" "cg11024682" "cg06500161" "cg27243685"

#print unique pubmed ids
unique(HDL_in_cat$PMID)
#[1] 28194238 25583993

#read in EWAS cat with total cholesterol studies
total_cat <- read.csv("Total_cholesterol_EWAS_cat.csv")

#check data dimensions 
dim(total_cat)
#[1] 250  33

#filter cat to whole blood, sample size > 1000 and P < 3.8e-8
total_cat <- total_cat[total_cat$Tissue == "Whole blood" & total_cat$P < 3.6e-8 & total_cat$N > 1000, ]

dim(total_cat)
#[1] 77 33

unique(total_cat$PMID)
#[1] 28213390 28173150


#check data dimensions 
dim(total_cat)
#[1] 77 33

#check overlap between our study and EWAS cat
total_in_cat <- total_cat[total_cat$CpG %in% cpgs,]

#check data dimensions 
dim(total_in_cat)
#[1]  5 33

#check num of unique cpgs
length(unique(total_in_cat$CpG))
#[1] 2

#print unique cpg names
unique(total_in_cat$CpG)
#[1] "cg17901584" "cg09737197"

#print unique pubmed ids 
unique(total_in_cat$PMID)
#[1] 28213390

#read in EWAS cat with WHR studies
WHR_cat <- read.csv("WHR_ewas_cat.csv")

#check data dimensions
dim(WHR_cat)
#[1] 10 33

#filter cat to whole blood, sample size > 1000 and P < 3.8e-8
WHR_cat <- WHR_cat[WHR_cat$Tissue == "Whole blood" & WHR_cat$P < 3.6e-8 & WHR_cat$N > 1000, ]

#check data dimensions
dim(WHR_cat)
#[1] 10 33

#print unique pubmed ids
unique(WHR_cat$PMID)
#[1] 32901515

#check overlap between our study and EWAS cat
WHR_in_cat <- WHR_cat[WHR_cat$CpG %in% cpgs,]

#check data dimensions 
dim(WHR_in_cat)
#[1]  3 33

#check num of unique cpgs
length(unique(WHR_in_cat$CpG))
#[1] 3

#print unique cpg names
unique(WHR_in_cat$CpG)
#[1] "cg19693031" "cg00574958" "cg06500161"

#print unique pubmed ids
unique(WHR_in_cat$PMID)
#[1] 32901515



#make lists of unique pubmed ids for all traits
w <-unique(WHR_in_cat$PMID)
t <- unique(total_in_cat$PMID)
h <- unique(HDL_in_cat$PMID)
g <- unique(Glucose_in_cat$PMID)
b <- unique(bmi_in_cat$PMID)


#combine lists
all <- c(w, t,  h, g, b)

#check unique pubmed ids for all traits to check for GS data 
unique(all)
#[1] 32901515 28213390 28194238 25583993 27019061 28095459 28002404 26110892
#[9] 25935004 26119815 31510868 29278407 29099282


### make table of cpgs to show where they overlap between traits
d1 <- as.data.frame(unique(bmi_in_cat$CpG))
dim(d1) #[1] 11  1

d1$Trait <- rep("BMI", 11)
names(d1)[1] <- "CpG"


d2 <- as.data.frame(unique(WHR_in_cat$CpG))
dim(d2) #[1] 3  1

d2$Trait <- rep("WHR", 3)
names(d2)[1] <- "CpG"


d3 <- as.data.frame(unique(HDL_in_cat$CpG))
dim(d3) #[1] 4  1

d3$Trait <- rep("HDL", 4)
names(d3)[1] <- "CpG"

d4 <- as.data.frame(unique(total_in_cat$CpG))
dim(d4) #[1] 2  1

d4$Trait <- rep("Total", 2)
names(d4)[1] <- "CpG"

d5 <- as.data.frame(unique(Glucose_in_cat$CpG))
dim(d5)

d5$Trait <- rep("Glucose", 2)
names(d5)[1] <- "CpG"

check <- rbind(d1,d2)
check <- rbind(check, d3)
check <- rbind(check, d4)
check <- rbind(check, d5)


table(check$CpG, check$Trait)

             BMI Glucose HDL Total WHR
cg00574958   1       0   0     0   1
cg06500161   1       1   1     0   1
cg08309687   1       0   0     0   0
cg09737197   0       0   0     1   0
cg10975897   1       0   0     0   0
cg11024682   1       1   1     0   0
cg14476101   1       0   0     0   0
cg16246545   1       0   0     0   0
cg17058475   1       0   0     0   0
cg17901584   1       0   1     1   0
cg19693031   0       0   0     0   1
cg23032421   1       0   0     0   0
cg27243685   1       0   1     0   0



