# This file contains functions to calculate gene likelihoods assuming an individual beta value between genes. 
# To run the calculation as a whole execute the function calculateLLIndivBeta. The total likelihood result
# will be returned from the function. The full results will be stored to an .RData file; to view these results 
# call load("./results/llindivbetaresults.RData"). 
# Author: Rori Rohlfs, Lars Gronvold, John Mendoza
# Date 2/25/19

#Include all our necessary libraries
library(mvtnorm)
source('./scripts/eve-io.R')

#' Estimate initial parameter values without using phylogeny
#'
#' @param gene.data A matrix of expression values with samples in columns and genes in rows
#' @param colSpecies A character vector with same length as columns in gene.data, specifying the 
#' species for the corresponding column.
#'
#' @return A matrix of initial parameters with the parameters in columns and genes in rows
initialParams <- function(gene.data, colSpecies)
{
  colSpeciesIndices <- split(seq_along(colSpecies), f = colSpecies)
  species.mean <- sapply(colSpeciesIndices, function(i){ rowMeans(gene.data[,i]) })
  species.var <- sapply(colSpeciesIndices, function(i){ apply(gene.data[,i],1,var) })
  
  theta <- rowMeans(species.mean)
  sigma2 <- apply(species.mean,1,var)
  alpha <- .5
  beta <- rowMeans(species.var) / sigma2
  
  return(cbind(theta,sigma2,alpha,beta))
}

# Expectation and variance from phylogeny
# Use this function to calculate expected species mean, and evolutionary
# variance for each node, given a common theta, alpha, and sigma^2 for all edges
calcExpVarOU <- function(tree, theta, alpha, sigma.squared, 
                         rootVar, rootE)
{
  # Declare vectors of expectation values and variances for all tips
  # and internal nodes
  expected.mean <- numeric( Ntip(tree) + Nnode(tree) )
  evol.variance <- numeric( Ntip(tree) + Nnode(tree) )
  
  # Get the order of the edges to traverse from root to tips
  edgeOrder <- rev(postorder(tree))
  
  # Get index of root
  root.index <- tree$edge[edgeOrder[1], 1]
  
  # Set root values
  expected.mean[root.index] = rootE
  evol.variance[root.index] = rootVar
  
  for(i in edgeOrder)
  {
    parent.index <- tree$edge[i, 1]
    child.index <- tree$edge[i, 2]
    delta.T <- tree$edge.length[i]
    expected.mean[child.index] <- expected.mean[parent.index] * exp(-alpha * delta.T) +
      theta * (1 - exp(-alpha * delta.T))
    evol.variance[child.index] <- evol.variance[parent.index] * exp(-2 * alpha * delta.T) +
      sigma.squared / (2 * alpha) * (1 - exp(-2 * alpha * delta.T))
  }
  
  return(data.frame(expected.mean, evol.variance))
}

# Covaraince from phylogeny

# Use this function to calculate the covariance matrix given alpha for each edge
calcCovMatOU <- function(tree, alpha, evol.variance)
{
  # Copy the phylogeny and multiply edge lengths with alpha
  attenuationTree <- tree
  attenuationTree$edge.length <- attenuationTree$edge.length * alpha
  
  # calculate the attenuation matrix using the cophenetic distance function
  attenuation.Matrix <- cophenetic(attenuationTree)
  
  # get matrix of variances of the most recent common ancestors
  variance.MRCA <- apply(mrca(tree), 1:2, function(i) evol.variance[i])
  
  # Return the covariance matrix
  return (variance.MRCA * exp(-attenuation.Matrix))
}

# Use this function to add within-species variance. This will expand
# the covariance matrix and species mean vector to account for within-species variance
expandECovMatrix <- function(expected.mean, covar.matrix, sigma.squared, alpha, beta, index.expand )
{
  # Expand covariance matrix
  covar.matrix.exp <- covar.matrix[index.expand, index.expand]
  # Add within species variance
  diag(covar.matrix.exp) <- diag(covar.matrix.exp) + beta * sigma.squared / (2 * alpha)
  
  # Expand expected values
  expected.mean.exp <- expected.mean[index.expand]
  
  return(list(cov.matr = covar.matrix.exp, expected.mean = expected.mean.exp))
}

# Maximum likelihood estimations

logLikOU <- function(param.matrix.row, tree, gene.data.row, index.expand)
{
  # Create varaibles with paramters to pass into the function for calculating expression variance 
  theta <- param.matrix.row[1]
  sigma.squared <- param.matrix.row[2]
  alpha <- param.matrix.row[3]
  beta <- param.matrix.row[4]
  
  # Define the expected species mean for the root and the evolutionary variance for the root
  expected.species.mean.root <- theta
  evol.var.root <- sigma.squared / (2 * alpha)
  
  
  expression.var <- calcExpVarOU(tree, theta, alpha, sigma.squared, rootVar = evol.var.root, rootE = expected.species.mean.root)
  covar.matrix <- calcCovMatOU(tree, alpha, evol.variance = expression.var$evol.variance)
  
  expanded.matrix <- expandECovMatrix(expected.mean = expression.var$expected.mean,covar.matrix,
                                      sigma.squared, alpha, beta, index.expand)
  
  # Get log likelihood from the multivariate normal distribution density
  ll <- dmvnorm(x = gene.data.row, mean = expanded.matrix$expected.mean, sigma = expanded.matrix$cov.matr, log = TRUE )
  return(ll)
}

calculateLLPerGene <- function(param.matrix.row, tree, gene.data.row, index.expand)
{
  ll <- -logLikOU(param.matrix.row, tree, gene.data.row, index.expand)
  return(ll)
}

# Use this function to calculate the total likelihood 
calculateTotalLL <- function(ll.pergene)
{
  llTotal <- abs(prod(ll.pergene))
  return(llTotal)
}

# Test for divergence and diversity between genes within species assuming an individual beta value for each gene
calculateLLIndivBeta <- function(tree, gene.data, colSpecies = colnames(gene.data))
{
  #Calculate the per gene parameter matrix based on the gene data
  initial.param.matrix <- initialParams(gene.data, colSpecies)
  
  # match the column species with the phylogeny tip labels
  index.expand <- match(colSpecies, tree$tip.label)
  
  # Create a vector to hold the values of the likelihood per gene and a list to hold the return values from the call to optim
  ll.pergene.IndivBeta <- vector(mode = "numeric", length = nrow(initial.param.matrix))
  max.params <- list()
  
  # For each gene, optimize the parameters and store the resulting likelihood in the likelihood vector
  for(row in 1:length(ll.pergene.IndivBeta))
  {
    max.params <- optim(initial.param.matrix[row, ], fn = calculateLLPerGene, gr = NULL, tree, gene.data[row, ], index.expand,
                     method = "L-BFGS-B", lower = c(.000000000000001, .000000000000001, .000000000000001, .000000000000001))
    ll.pergene.IndivBeta[row] <- as.numeric(max.params[2])
  }
  
  # Calculate the total likelihood as the product of all values within the likelihood vector
  ll.total.IndivBeta <- calculateTotalLL(ll.pergene.IndivBeta)
  
  # Save the results to a file that can be loaded into any future R environment for future use
  # To load this data simply call load("./results/llindivbetaresults.RData")
  save(ll.pergene.IndivBeta, ll.total.IndivBeta, file = "./results/llindivbetaresults.RData")
  
  return(ll.total.IndivBeta)
}


