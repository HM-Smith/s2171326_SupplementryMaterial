library("biglasso")
library("bigmemory")
library("caret")
library("data.table")
library("dplyr")

#Load methylation data
load("20180128_0046_betas_ct_batch.Rdata")

#check data 
head(betas_ct_batch)[1:5, 1:5] # samples are columns and CpGs are rows

#sample names should be rows and CpGs columns so data needs transposing
meth_data = t(betas_ct_batch)

#check data dimensions
dim(meth_data)
#[1]   1469 255356

#Load outcome
response = read.csv("training_data_80.csv")


# Align IDs between dataframes 
meth_data=meth_data[which(row.names(meth_data) %in% response$Sample),]


probes = read.csv("probe_04_reliability.csv")


probes = probes$CpG

#subset to reliable probes
meth_data = meth_data[, colnames(meth_data) %in% probes]


# Align IDs between dataframes 
response=response[which(response$Sample %in% rownames(meth_data)),]


#make id list
ids=response$Sample

#match meth rownames to response
meth_data=meth_data[match(ids,row.names(meth_data)),]

identical(rownames(meth_data), response$Sample) # TRUE

# Convert data to big.matrix format (required for biglasso)
X <- as.big.matrix(as.matrix(meth_data))
y <- as.numeric(response$protein)  # Ensure response is numeric or binary

## function to put Twin IDs in same fold, this also means longitudinal measures
## will remain in same fold
getGroupCVFoldIDs <- function(groups, nFolds, seed = NULL) {
  if (!is.null(seed)) {
    set.seed(seed)
  }
  
  counts <- table(groups)
  elementsPerFold <- length(groups) / nFolds
  
  usedGroups <- character()
  currentFoldNumber <- 1
  groupToFoldMap <- list()
  foldToNElementsMap <- list()
  currentFoldNElements <- 0
  
  for (i in seq(length(counts))) {
    remainingGroups <- setdiff(names(counts), usedGroups)
    sampledGroup <- sample(remainingGroups, size = 1)
    usedGroups <- c(usedGroups, sampledGroup)
    currentFoldNElements <- currentFoldNElements + counts[[sampledGroup]]
    groupToFoldMap[[sampledGroup]] <- currentFoldNumber
    if ((currentFoldNElements >= elementsPerFold) || (i == length(counts))) {
      foldToNElementsMap[[currentFoldNumber]] <- currentFoldNElements
      if (currentFoldNumber < nFolds) {
        currentFoldNumber <- currentFoldNumber + 1
        currentFoldNElements <- 0
      }
    }
  }
  foldIDs <- sapply(groups, function(x) {
    groupToFoldMap[[as.character(x)]]
  })
  
  list(foldIDs = foldIDs, 
       groupToFoldMap = groupToFoldMap, 
       foldToNElementsMap = foldToNElementsMap)
}

#apply fold group function to our data 
folds <- getGroupCVFoldIDs(response$PAIRID, 5, 123)

#check spread across folds
table(folds$foldIDs)


#make list of fold assignments 
fold_ids <- folds$foldIDs

# Run cross-validation with grouping
lasso.cv <- cv.biglasso(
  X, y,
  family = "gaussian",
  alpha = 0.5,
  ncores = 2,
  nfolds = 5,
  cv.ind = fold_ids  # Pass fold assignments to cv.biglasso
)

# Fit the model using the best lambda from cross-validation
fit <- biglasso(
  X, y,
  family = "gaussian",
  alpha = 0.5,
  ncores = 2,
  lambda = lasso.cv$lambda.min
)


#extract betas
coefs <- coef(fit)
#set as dataframe 
coefs1=as.data.frame(as.matrix(coefs))
#make cpg column
coefs1$CpG=row.names(coefs1)
#rename beta column
names(coefs1)[1]="Coefficient"
#subset to non-zero coefficients 
coefs2=coefs1[which(coefs1$Coefficient!=0),]
#remove 
coefs2=coefs2[-1,]
#create an absolute coefficient column
coefs2$Rank=abs(coefs2$Coefficient)
#check number of CpGs
dim(coefs2)



# Save coefficients
saveRDS(coefs2, file = "") 
