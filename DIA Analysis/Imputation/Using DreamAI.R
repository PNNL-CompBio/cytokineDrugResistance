require("cluster")
require("survival")
require("randomForest")
require("missForest")
require("glmnet")
require("Rcpp")
require("foreach")
require("itertools")
require("iterators")
require("Matrix")
require("impute")
require("DreamAI")

library(kableExtra)
library(gridExtra)
library(dplyr)
library(MSnSet.utils)
library(ggplot2)
library(clusterProfiler)
library(amlresistancenetworks)
source("../util/synapseUtil.R")
source("../util/loadData.R")

## Using synapse ID to fecth table + metadata
dat <- load.DIA.data("syn26164985")
prot.dat <- dat[["Long form data"]]
meta <- dat[["Metadata"]] 
colnames(meta)[4] <- "Trametinib"
rownames(meta) <- meta$Sample

## Turning long form into matrix
## for any peptides with the same parent protein we average intensity.
prot.mat <- prot.dat %>%
  select(Protein, Sample, Intensity) %>%
  tidyr::pivot_wider(values_from='Intensity',names_from='Sample',
                     values_fn=list(Intensity = mean)) %>%
  tibble::column_to_rownames('Protein') %>%
  as.matrix()

data <- prot.mat
use_MissForest <- T

# remove rows where all data is NA
data = data[rowSums(is.na(data)) != ncol(data), ]

# Any features with more than 90% missing data are thrown out. Needed for using
# DreamAI!! using a slightly stricter filter in order to obtain under 90% missing
# in the columns as well!!
N <- ncol(data)
keeping <- rowSums(is.na(data)) < 0.9*N 
data <- data[keeping, ]


0.9*nrow(data)
colSums(is.na(data))
dropped <- 25
i = 0

while(any(colSums(is.na(data)) >= 0.9*nrow(data))){
  weights <- colSums(is.na(data))/(0.9*nrow(data))
  weights <- matrix(weights, nrow = nrow(data), ncol = 30, byrow = TRUE)
  
  score <- as.numeric(!is.na(data))*weights
  rownames(score) <- rownames(data)
  importance <- rowSums(score) %>%
    sort(decreasing = F) %>%
    names()
  
  keeping <- importance[dropped:length(importance)]
  data <- data[keeping, ]
  print(i)
  i <- i + 1
}


if (use_MissForest) {
  impute<- DreamAI(data,k=10,maxiter_MF = 10, ntree = 100,maxnodes = NULL,maxiter_ADMIN=30,tol=10^(-2),gamma_ADMIN=0,gamma=50,CV=FALSE,fillmethod="row_mean",maxiter_RegImpute=10,conv_nrmse = 1e-6,iter_SpectroFM=40, method = c("KNN", "ADMIN", "MissForest", "Brinn", "SpectroFM", "RegImpute"),out="Ensemble")
} else {
  impute<- DreamAI(data,k=10,maxiter_MF = 10, ntree = 100,maxnodes = NULL,maxiter_ADMIN=30,tol=10^(-2),gamma_ADMIN=0,gamma=50,CV=FALSE,fillmethod="row_mean",maxiter_RegImpute=10,conv_nrmse = 1e-6,iter_SpectroFM=40, method = c("KNN", "ADMIN", "Brinn", "SpectroFM", "RegImpute"),out="Ensemble")
}

write.table(impute$Ensemble,"imputed_file.tsv",row.names=TRUE, col.names=NA, sep='\t', quote = FALSE)

corre <- sapply(rownames(xx), function(name){
  cor(xx[name,], zz[name,])
})


