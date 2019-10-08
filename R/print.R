#' @title Print method for pcpr2
#' @description
#' \code{print.pcpr2} Prints a mimimal output.
#' @details
#' This function prints the model fit and Partial R2 values calculated
#' by the \code{\link{runPCPR2}} function.
#' @param x an object of class \code{pcpr2}.
#' @param ... additional arguments passed to the function.
#' @method print pcpr2
#' @export
#' @return the input object is returned silently.
#' @examples
#' output <- runPCPR2(transcripts, Z_metadata)
#' print(output)
print.pcpr2 <-
  function(x, ...){
    if (!inherits(x, "pcpr2"))
      stop("Object must be of class 'pcpr2'")
    cat("Principal Component R-squared method (PC-PR2) \n")
    cat("Number of observations:", x$dimensions[1], "Number of X-variables:", x$dimensions[2], "\n")
    cat("Partial R2 for Z-variables:\n")
    print(x$pR2, ...)
  }
