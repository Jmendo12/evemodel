


# test for shift in theta at branch specified by the extant species "shiftSpecies"
twoThetaTest <- function(tree, gene.data, isTheta2edge, colSpecies = colnames(gene.data), ...){

  cat("Fit with single theta...\n") 
  oneThetaRes <- fitOneTheta(tree, gene.data, ...)
  cat("Fit with two thetas...\n")
  twoThetaRes <- fitTwoTheta(tree, gene.data,isTheta2edge = isTheta2edge, ...)
  
  LRT <- 2 * (twoThetaRes$ll - oneThetaRes$ll)
  
  return( list(oneThetaRes = oneThetaRes,
               twoThetaRes = twoThetaRes,
               LRT = LRT) )
}
