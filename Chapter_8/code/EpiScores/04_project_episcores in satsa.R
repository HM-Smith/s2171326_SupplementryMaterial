# Load requisite libraries 
library(data.table)

#Load methylation data
load("20180128_0046_betas_ct_batch.Rdata")


#check dimensions 
dim(betas_ct_batch)


#transpose 
meth = t(betas_ct_batch)

#set as dataframe 
meth = as.data.frame(meth)

### read in weights ###   
g_wts <- readRDS("")

#check dimensions 
dim(g_wts)

#select cpg and betas 
meanBetas = g_wts[,c(2,1)]
names(meanBetas)=c("CpG","Beta")


### subset to relevant CpG sites ###
a36 = which(names(meth) %in% meanBetas$CpG)

#subset meth to cpgs in episcore
meth <- meth[,a36]

#check for na
table(is.na(meth))
# 
# 
# FALSE
# 20566


### line up weights and CpGs ###
b36 = which(meanBetas$CpG %in% colnames(meth)) 
mean_betas <- meanBetas[b36,]
meth <- meth[,match(colnames(meth), mean_betas$CpG)]
meth <- meth[,mean_betas$CpG]

identical(names(meth), mean_betas$CpG) 
#[1] TRUE
##covert methylation to matrix
meth = as.matrix(meth)

### create predictor ###
pred36 <- meth %*% mean_betas$Beta
pred36 <- as.data.frame(pred36)
names(pred36) <- "X_EpiScore"
pred36$Sample <- rownames(pred36)

#check dimensiosn 
dim(pred36)

#save out episcores
saveRDS(pred36, file = "")
