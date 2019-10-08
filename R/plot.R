#' @title Plot PC-PR2 output
#' @description
#' \code{plot.pcpr2} plots partial R2 values.
#' @details
#' This function plots a barchart of partial R2 values
#' created by the \code{\link{pcpr2}} function for easy interpretation.
#' The bars are annotated with values.
#' @param x an object of class \code{pcpr2}.
#' @param ... additional arguments passed to the
#' \code{\link{barplot}} function.
#' @method plot pcpr2
#' @export
#' @return NULL
#' @examples
#' output <- runPCPR2(transcripts, Z_metadata)
#' plot(output, col="red", main="Variability in transcriptomics data explained by covariates",
#' ylab="Rpartial2")
plot.pcpr2 <- function(x, ...){
  if (!inherits(x, "pcpr2"))
    stop("Object must be of class 'pcpr2'")
  data <- x$pR2
  opar <- par(no.readonly=TRUE)
  par(mar=c(6,5,4,2))
  bp <- barplot(unname(data),
          xlab = "",
          ylab = "Weighted Rpartial2",
          ylim = c(0, max(data) * 1.3),
          las = 2, ...)
  axis(1, at = bp, labels = c(names(data)), cex.axis = 0.8, las=2)
  rounded <- round(data, 3)
  text(bp, data, labels = rounded, pos = 3, cex = 0.8)
  par(opar)
}
