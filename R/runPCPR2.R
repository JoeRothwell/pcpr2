#' The Principal Component R-squared method (PC-PR2)
#'
#' Runs the PC-PR2 method on X, an omics matrix of intensities, and Z, the subject metadata to be assessed in the model.
#' @param X A matrix of omics data. It is recommended that the data be appropriately scaled and log transformed if necessary.
#' @param Z A data frame of subject metadata, with the same number of rows as X, whose influence on the omics data is to be assessed. Categorical variables should be coded as factors.
#' @param pct.threshold The proportion of variability desired to be explained. Defaults to 0.8.
#' @keywords pcpr2, Rpartial2, principal component analysis, pca, omics, metabolomics, transcriptomics
#' @examples output <- runPCPR2(transcripts, Z_metadata)
#' @examples output$pR2
#' @export
#' @return Returns a list with 6 elements:
#' \item{dimensions}{dimensions of the -omics matrix}
#' \item{types}{data types for each Z-variable}
#' \item{threshold}{the chosen variability to be explained}
#' \item{pcn}{the number of principal components used for the chosen threshold}
#' \item{anovaFit}{ANOVA output for each PC}
#' \item{pR2}{Rpartial2 values for each Z-variable and the overall R2}
runPCPR2 <- function(X, Z, pct.threshold = 0.8) {

  # The following commented lines are not used but have been kept from the original script.
  # Scaling and centering should now be done outside the analysis
  # Load the data in the X matrix containing NMR spectra and Z metadata containing the list of explanatory variables of interest
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
  X_DataMatrix <- X
  Z_Meta <- Z

  # Check arguments:
  # X is a numeric matrix with no NAs
  stopifnot(is.numeric(X_DataMatrix),
            is.matrix(X_DataMatrix),
            !anyNA(X_DataMatrix),
  # Same number of n for X and metadata
            nrow(X_DataMatrix) == nrow(Z_Meta))

  if(ncol(X_DataMatrix) < 4) {
    warning(paste("Only", ncol(X_DataMatrix),
                  "column(s) in your X matrix. Are you sure you need the PC-PR2?"))
  }

  Z_MetaRowN <- nrow(Z_Meta)
  Z_MetaColN <- ncol(Z_Meta)
  ColNames   <- names(Z_Meta)
  #ColNames1  <- c(ColNames, "Rmodel2")

  # Obtain eigenvectors (the following hashed code not efficient for long matrices)
  #X_DataMatrix_t <- t(X_DataMatrix)
  #symMat <- X_DataMatrix %*% X_DataMatrix_t
  #eigenData    <- eigen(symMat)
  #eigenVal     <- eigenData$values
  #eigenVecMat  <- eigenData$vectors
  #percents_PCs <- eigenVal/sum(eigenVal)

  # Set variability to be explained and get number of PCs required (force min to 3)
  #my_counter_2 <- sum(1 - cumsum(rev(percents_PCs)) <= pct.threshold)
  #if(my_counter_2 > 3) pc_n <- my_counter_2 else pc_n <- 3

  #pc_data_matrix <- eigenVecMat[, 1:pc_n ]

  # The following (contributed by Vivian) is much faster for long matrices
  respca         <- prcomp(X_DataMatrix, center = F, scale. = F)
  pc_n           <- max(3, min(which(cumsum(respca$sdev^2)/ncol(X_DataMatrix) >= pct.threshold)))
  eigenVal       <- respca$sdev^2
  pc_data_matrix <- respca$x[, 1:pc_n ]

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
  AnovaTab <- car::Anova(multifit, type=3, singular.ok = F)
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
  names(pR2Sums) <- c(ColNames, "R2")

  # Return results
  output <- list(dimensions = dim(X_DataMatrix),
                 types = sapply(Z_Meta, class),
                 threshold = pct.threshold,
                 pcn = pc_n,
                 anovaFit = multifit,
                 pR2 = pR2Sums)
  class(output) <- c("pcpr2", "list")
  return(output)
}

