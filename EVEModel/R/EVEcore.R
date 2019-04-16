# This file contains the core functions to calculate the likelihood of the observed gene expression
# given a set of parameters for the EVE model
# Author: Rori Rohlfs, Lars Gronvold, John Mendoza
# Date 2/25/19

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

#' Expectation and variance from phylogeny under the OU model
#'
#' Calculate expected species mean, and evolutionary variance for each node in the tree
#' given theta, sigma^2 and alpha parameters for each edge in the tree.
#'
#' @param tree Species phylogeny (phylo object)
#' @param thetas Vector of theta parameters for each edge in the tree
#' @param alphas Vector of alpha parameters for each edge in the tree
#' @param sigma2s Vector of sigma2 parameters for each edge in the tree
#' @param rootVar The evolutionary variance at the root node
#' @param rootE The expected species mean at the root node
#'
#' @return data.frame with expected.mean and evol.variance for each node
#' @export
#'
#' @examples
calcExpVarOU <- function(tree, thetas, alphas, sigma2s,
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
    expected.mean[child.index] <- expected.mean[parent.index] * exp(-alphas[i] * delta.T) +
      thetas[i] * (1 - exp(-alphas[i] * delta.T))
    evol.variance[child.index] <- evol.variance[parent.index] * exp(-2 * alphas[i] * delta.T) +
      sigma2s[i] / (2 * alphas[i]) * (1 - exp(-2 * alphas[i] * delta.T))
  }

  return(data.frame(expected.mean, evol.variance))
}

#' Covariance from phylogeny under the OU model
#'
#' Calculates the OU model covariance matrix given alpha parameters for each edge
#' and the evolutionary variance at each node of a tree.
#'
#' @param tree Species phylogeny (phylo object)
#' @param alphas Vector of alpha parameters for each edge in the tree
#' @param evol.variance
#'
#' @return Covariance matrix
#' @export
#'
#' @examples
calcCovMatOU <- function(tree, alphas, evol.variance)
{
  # Copy the phylogeny and multiply edge lengths with alpha
  attenuationTree <- tree
  attenuationTree$edge.length <- attenuationTree$edge.length * alphas

  # calculate the attenuation matrix using the cophenetic distance function
  attenuation.Matrix <- cophenetic(attenuationTree)

  # get matrix of variances of the most recent common ancestors
  variance.MRCA <- apply(mrca(tree), 1:2, function(i) evol.variance[i])

  # Return the covariance matrix
  return (variance.MRCA * exp(-attenuation.Matrix))
}

# Use this function to add within-species variance. This will expand
# the covariance matrix and species mean vector to account for within-species variance
expandECovMatrix <- function(expected.mean, covar.matrix, sigma2, alpha, beta, index.expand )
{
  # Expand covariance matrix
  covar.matrix.exp <- covar.matrix[index.expand, index.expand]
  # Add within species variance
  diag(covar.matrix.exp) <- diag(covar.matrix.exp) + beta * sigma2 / (2 * alpha)

  # Expand expected values
  expected.mean.exp <- expected.mean[index.expand]

  return(list(cov.matr = covar.matrix.exp, expected.mean = expected.mean.exp))
}

# silly function with same values of parameters across the tree
calcExpVarOUconst <- function(tree, theta, alpha, sigma2){
  N <- Nedge(tree)
  calcExpVarOU(tree, thetas = rep(theta,N), alphas = rep(alpha,N), sigma2s = rep(sigma2,N),
               rootVar = sigma2/(2*alpha), rootE = theta)
}

# Likelihood of gene expression with given parameters
logLikOU <- function(theta, sigma2, alpha, beta, tree, gene.data.row, index.expand)
{
  # Define the expected species mean for the root and the evolutionary variance for the root
  expected.species.mean.root <- theta
  evol.var.root <- sigma2 / (2 * alpha)


  expression.var <- calcExpVarOUconst(tree = tree, theta = theta, alpha = alpha, sigma2 = sigma2)
  covar.matrix <- calcCovMatOU(tree, alphas = alpha, evol.variance = expression.var$evol.variance)

  expanded.matrix <- expandECovMatrix(expected.mean = expression.var$expected.mean,covar.matrix,
                                      sigma2, alpha, beta, index.expand)

  # Get log likelihood from the multivariate normal distribution density
  ll <- dmvnorm(x = gene.data.row, mean = expanded.matrix$expected.mean, sigma = expanded.matrix$cov.matr, log = TRUE )
  return(ll)
}

# Likelihood of gene expression with given parameters
# isThetaShiftEdge logical vector specifying if the theta2 or theta1 shall be used
logLikTwoTheta <- function(theta1, theta2, sigma2, alpha, beta, tree, isThetaShiftEdge, gene.data.row, index.expand)
{
  N <- Nedge(tree)
  expression.var <- calcExpVarOU(tree, thetas = ifelse(isThetaShiftEdge,theta2,theta1),
                                 alphas = rep(alpha,N), sigma2s = rep(sigma2,N),
                                 rootVar = sigma2/(2*alpha), rootE = theta1)

  covar.matrix <- calcCovMatOU(tree, alphas = alpha, evol.variance = expression.var$evol.variance)

  expanded.matrix <- expandECovMatrix(expected.mean = expression.var$expected.mean,covar.matrix,
                                      sigma2, alpha, beta, index.expand)

  # Get log likelihood from the multivariate normal distribution density
  ll <- dmvnorm(x = gene.data.row, mean = expanded.matrix$expected.mean, sigma = expanded.matrix$cov.matr, log = TRUE )
  return(ll)
}
