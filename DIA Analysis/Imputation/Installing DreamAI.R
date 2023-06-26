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
require("devtools")
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("impute")
require("impute")


## Converting warning to error, so use this?
Sys.setenv(R_REMOTES_NO_ERRORS_FROM_WARNINGS= "true")

require("remotes")
install_github("WangLab-MSSM/DreamAI/Code")
