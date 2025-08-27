# Load packages and functions for use in models 
library(bigmemory)
library(optparse)
library(biglasso)
library(data.table)
library(foreign)

# Read in Methylation dataframe, transponse and add in cpg names 
methx=fread("gs20k-betas-scaled-cleaned.csv")
methx=as.data.frame(methx)
methx=t(methx)
methx=as.data.frame(methx)
cpgs=readRDS("CpGs_combined.rds") # read in cpg names
names(methx)=cpgs

# Read in phenotype for prediction
phen=read.table("whr_res.phen",header=T)
# Rename phenotype column for convenience 
names(phen)[3]="phen"
# Remove NAs for this analysis 
phen=phen[which(!is.na(phen$phen)),]

# Align IDs between dataframes 
methx=methx[which(row.names(methx) %in% phen$FID),]
ids=phen$FID
methx=methx[match(ids,row.names(methx)),]


################################################################################

# Split big lasso matricies using cbindBM - load relevant functions for this

################################################################################

# x needs to be a list of big matrices, and they have to have the same rownames
cbindBM_list <- function(x, binding="right", 
                         z=NULL, type=NULL, separated=NULL,
                         backingfile=NULL, backingpath=NULL,
                         descriptorfile=NULL, binarydescriptor=FALSE,
                         shared=TRUE, erase = TRUE)
{
  
  if (is.null(type)) type <- typeof(x[[1]])
  if (is.big.matrix(x[[1]])) {
    if (is.null(separated)) separated <- is.separated(x[[1]])
  } else {
    separated <- FALSE
  }
  
  cols_list <- list()
  total_cols <- 0
  for (i in 1:length(x)) {
    cols <- cleanupcols(NULL, ncol(x[[i]]), colnames(x[[i]]))
    cols_list <- append(cols_list, list(cols))
    total_cols <- total_cols + ncol(x[[i]])
  }    
  
  if (is.null(z)) {
    z <- big.matrix(nrow=nrow(x[[1]]), ncol=total_cols, type=type, init=NULL,
                    dimnames=dimnames(x[[1]]), separated=separated,
                    backingfile=backingfile, backingpath=backingpath,
                    descriptorfile=descriptorfile,
                    binarydescriptor=binarydescriptor, shared=shared)
  }
  
  counter <- 0
  for (i in 1:length(cols_list)) {
    print(i)
    if (i == 1) {
      z[, 1:length(cols_list[[i]])] <- x[[i]][,cols_list[[i]]]
    } else {
      z[, (counter + 1):(counter + length(cols_list[[i]]))] <- x[[i]][,cols_list[[i]]]
    }
    counter <- counter + length(cols_list[[i]])
    print(counter)
    
    if (erase == TRUE) {
      cat("\nErasing chunk and liberating memory...\n\n")
      x[[i]] <- "Replacement"
      gc()
    }
  }
  return(z)
}


cleanupcols <- function(cols=NULL, nc=NULL, colnames=NULL) {
  if (is.null(cols)) cols <- 1:nc
  else {
    if (!is.numeric(cols) & !is.character(cols) & !is.logical(cols))
      stop("column indices must be numeric, logical, or character vectors.")
    if (is.character(cols))
      if (is.null(colnames)) stop("column names do not exist.")
    else cols <- mmap(cols, colnames)
    if (is.logical(cols)) {
      if (length(cols) != nc)
        stop(paste("column vector length must match the number of",
                   "columns of the matrix."))
      cols <- which(cols)
    }
    tempj <- .Call("CCleanIndices", as.double(cols), as.double(nc), PACKAGE="bigmemory")
    if (is.null(tempj[[1]])) stop("Illegal column index usage in extraction.\n")
    if (tempj[[1]]) cols <- tempj[[2]]
  }
  return(cols)
}

meth=methx
div <- 30 # Number of chunks to divide OG methylation dataframe

por <- ceiling(length(colnames(meth))/div)
chunk_list <- list()

for (i in 1:div) {
  cat(paste0("\nWorking on chunk: ", i, " of ", div))
  if (i == 1) {
    chunk <- as.big.matrix(meth[,1:(por-1)])
  } else if (i == div) {
    chunk <- as.big.matrix(meth[,(por*(i-1)):length(colnames(meth))])
  } else {
    chunk <- as.big.matrix(meth[,(por*(i-1)):((por*i)-1)])
  }
  cat("\nMade chunk. Appending to chunk list...\n")
  chunk_list <- append(chunk_list, list(chunk))
  gc()
}

# Saving names prior to chunk fusing
names <- colnames(meth)
rm(meth)

cat("\nRAM clean up...\n\n")
gc()

cat("\nFusing chunks!\n\n")
methx <- cbindBM_list(x = chunk_list)
rm(chunk, chunk_list)

# Set CpG names
options(bigmemory.allow.dimnames=TRUE)
colnames(methx)<- names



set.seed(1783) # set seed to ensure fold variation minimised 

identical(rownames(methx), phen$FID) # TRUE

# Assign ytrain
y <- as.numeric(phen$phen) # Create a numeric list for this variable to feed in as y to the model

# Run training 
lasso.cv <- cv.biglasso(methx, y, family="gaussian", alpha = 0.5, ncores = 8, nfolds = 20) # cross validation to get best lambda
fit <- biglasso(methx, y, family = "gaussian", alpha = 0.5, ncores = 8, lambda = lasso.cv$lambda.min) # model fit 
coefs <- coef(fit) # Get betas
coefs1=as.data.frame(as.matrix(coefs))
coefs1$CpG=row.names(coefs1)
names(coefs1)[1]="Coefficient"
coefs2=coefs1[which(coefs1$Coefficient!=0),]
coefs2=coefs2[-1,]
coefs2$Rank=abs(coefs2$Coefficient)
# Save coefficients
saveRDS(coefs2, "WHR_predictor.rds") 
