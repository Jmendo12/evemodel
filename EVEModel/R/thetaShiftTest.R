# This file contains the functions for calculating likelihoods using the theta shift test.
# Authors: Rori Rohlfs, Lars Gronvold, John Mendoza
# Date: 3/19/19

thetaShiftTest <- function(tree, geneData, shiftedSpecies, colSpecies = colnames(geneData))
{
  oneThetaRes <- fitIndivBeta(tree, geneData, colSpecies)

  twoThetaParams <- initialParamsTwoThetas(geneData, colSpecies, shiftedSpecies)

}

initialParamsTwoThetas <- function(geneData, colSpecies, shiftedSpecies)
{
  colSpeciesIndices <- split(seq_along(colSpecies), f = colSpecies)
  speciesMean <- sapply(colSpeciesIndices, function(i) { rowMeans(geneData[, i]) })
  speciesVar <- sapply(colSpeciesIndices, function(i) { apply(geneData[,i], 1, var) })

  theta1 <- rowMeans(speciesMean[, !(colnames(speciesMean) %in% shiftedSpecies)])
  theta2 <- rowMeans(speciesMean[, shiftedSpecies])
  sigma2 <- apply(speciesMean, 1, var)
  alpha <- .5
  beta <- rowMeans(speciesVar) / sigma2

  return(cbind(theta1, theta2, sigma2, alpha, beta))
}
