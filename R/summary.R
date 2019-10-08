#' @title Summary method for pcpr2
#' @description
#' \code{print.pcpr2} Prints a detailed summary of the PC-PR2 object
#' @details
#' This summary function prints a formatted summary the data supplied, the parameters supplied,
#' the model fit and the Partial R2 values calculated by the \code{\link{runPCPR2}} function.
#' @param x an object of class \code{pcpr2}.
#' @param ... additional arguments passed to the function.
#' @method summary pcpr2
#' @export
#' @return the input object is returned silently.
#' @examples
#' output <- runPCPR2(transcripts, Z_metadata)
#' summary(output)
summary.pcpr2 <- function(x, ...){
    if (!inherits(x, "pcpr2"))
      stop("Object must be of class 'pcpr2'")
    cat("Principal Component R-squared method (PC-PR2) summary \n\n")
    cat("X-matrix dimensions:", x$dimensions, "\n")
    cat("Z-variables, n =", length(x$types), ":\n")
    print(data.frame(class = x$types))
    cat("\n")
    cat("Proportion of variability desired to be explained:", x$threshold, "\n")
    cat("Principal components required:", x$pcn, "\n\n")
    cat("Overall and partial R2 for Z-variables:\n")
    print(data.frame(percentage = x$pR2), ...)
  }
