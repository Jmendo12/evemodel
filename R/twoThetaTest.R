


#' Test for shift in theta at branch
#'
#' @param tree Species phylogeny (phylo object)
#' @param gene.data A matrix of expression values with samples in columns and genes in rows
#' @param isTheta2edge Logical vector with same length as number of edges in tree specifying whether
#' the corresponding edge theta parameter should be theta2 (TRUE) or theta1 (FALSE)
#' @param colSpecies A character vector with same length as columns in gene.data, specifying the
#' species for the corresponding column.
#' @param ... parameters passed to \code{\link{fitOneTheta}} and \code{\link{fitTwoTheta}}
#'
#' @return List with: indivBetaRes (results from \code{\link{fitOneTheta}})
#' \itemize{
#'   \item oneThetaRes: results from \code{\link{fitOneTheta}}
#'   \item twoThetaRes: results from \code{\link{fitTwoTheta}}
#'   \item LRT: log likelihood ratio test statistic between the two-theta and one-theta model
#' }
#' @export
twoThetaTest <- function(tree, gene.data, isTheta2edge, colSpecies = colnames(gene.data), ...){

  cat("Fit with single theta...\n") 
  oneThetaRes <- fitOneTheta(tree, gene.data, colSpecies = colSpecies, ...)
  cat("Fit with two thetas...\n")
  twoThetaRes <- fitTwoTheta(tree, gene.data,isTheta2edge = isTheta2edge, colSpecies = colSpecies, ...)
  
  LRT <- 2 * (twoThetaRes$ll - oneThetaRes$ll)
  
  return( list(oneThetaRes = oneThetaRes,
               twoThetaRes = twoThetaRes,
               LRT = LRT) )
}
