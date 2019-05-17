#' The Principal Component R-squared method (PC-PR2)
#'
#' Runs PC-PR2 on X, an omics matrix of intensities, and Y, the subject metadata to be assessed in the model.
#' @param X A matrix of omics data. It is recommended that the data be approriately scaled and log transformed if necessary.
#' @param Y The metadata variables whose influence on matrix X is to be assessed. Variables should be converted to factors where necessary.
#' @param pct_threshold The proportion of variability desired to be explained. Defaults to 0.8.
#' @keywords pcpr2
#' @export
#' @examples runPCPR2(Xtestdata, Ytestdata)

runPCPR2 <- function(X, Y, pct_threshold = 0.8) {

  library(car)
  library(MetabolAnalyze)
  library(tidyverse)

  # Load the data in the X matrix containing NMR spectra and Z matrix containing the list of explanatory variables of interest
  #myPath <- "Documents/PCPR2/" Metabolomics_data <- "X_MetaboMatrix.TXT"
  #InterestFactors_data <- "Z_FactorMatrix.TXT"
  #Metabo_FilePath = paste(myPath,Metabolomics_data, sep="")
  #Factors_FilePath = paste(myPath,InterestFactors_data, sep="")
  #X_DataMatrix <- read.delim(Metabo_FilePath, row.names = 1, header = TRUE, sep = "\t")
  #Z_InterestFactors <- read.delim(Factors_FilePath, sep = "\t", header = TRUE, row.names = 1)

  # Center the data / Scale the data, edit the parameter "pareto" or "unit" of scaling according to your need
  #X_DataMatrixCentered = scale(X_DataMatrix, center = TRUE, scale = FALSE)
  #X_DataMatrixScaled = scaling(X_DataMatrixCentered, type = "pareto")

  # Load test data
  X_DataMatrixScaled <- X
  Z_Meta <- Y

  Z_MetaRowN <- nrow(Z_Meta)
  Z_MetaColN <- ncol(Z_Meta)
  ColNames   <- names(Z_Meta)

  # Obtain eigenvectors
  pct_threshold <- 0.8 # set variability desired to be explained
  X_DataMatrixScaled_t <- t(X_DataMatrixScaled)
  symMat <- X_DataMatrixScaled %*% X_DataMatrixScaled_t
  eigenData    <- eigen(symMat)
  eigenVal     <- eigenData$values
  eigenVecMat  <- eigenData$vectors
  percents_PCs <- eigenVal/sum(eigenVal)

  # Get number of PCs required for threshold (force min to 3)
  my_counter_2 <- sum(1 - cumsum(rev(percents_PCs)) <= 0.8)
  if(my_counter_2 > 3) pc_n <- my_counter_2 else pc_n <- 3

  pc_data_matrix <- eigenVecMat[, 1:pc_n ]

  #Perform linear multiple regression models on each eigenvector with factors of interest as explanatory variables
  #Categorical variables should be processed by as.factor, whereas continuous variables should not.
  #To be edited with your factors names

  # Convert categorical variables to factors (put them in varlist)
  #varlist <- c("sex", "smoking.status")
  #Z_Meta <- Z_Meta %>% mutate_at(vars(varlist), as.factor)

  DataCol <- Z_MetaColN + 1

  # Run a linear model with each eigenvector as the response
  TotSumSq <- apply(pc_data_matrix, 2, var) * (Z_MetaRowN - 1)
  multifit <- lm(pc_data_matrix ~ ., data = Z_Meta)

  # Run type 3 ANOVA on each PC
  AnovaTab <- Anova(multifit, type=3, singular.ok = F)
  SSP      <- AnovaTab$SSP

  # Extract sum of squares for each factor, removing intercept column
  # Need to take the diagonal of each factor matrix to get sums of squares
  Residuals  <- diag(AnovaTab$SSPE)
  RR         <- Residuals/TotSumSq

  type3mat0 <- sapply(SSP, diag)[, -1]
  type3mat  <- cbind(type3mat0, "SumSqResiduals" = Residuals)
  ST_ResidualR2 <- cbind("ST_R2" = 1-RR, "ST_Residuals" = RR)

  #Create partial R2 matrix and weighted matrix
  partialR2mat <- type3mat[, -DataCol] / (type3mat[, -DataCol] + type3mat[, DataCol])
  eigenVal     <- eigenVal[1:pc_n]
  weight       <- eigenVal/sum(eigenVal)

  partialR2MatWtProp <- cbind(partialR2mat, ST_ResidualR2[, 1]) * weight
  colnames(partialR2MatWtProp) <- NULL
  pR2Sums <- colSums(partialR2MatWtProp) * 100
}
