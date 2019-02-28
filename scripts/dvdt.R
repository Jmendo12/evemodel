# This file contains functions to calculate gene likelihoods assuming an individual beta value between genes. 
# To run the calculation as a whole execute the function calculateLLIndivBeta. The total likelihood result
# will be returned from the function. The full results will be stored to an .RData file; to view these results 
# call load("./results/llindivbetaresults.RData"). 
# Author: Rori Rohlfs, Lars Gronvold, John Mendoza
# Date 2/25/19

#Include all our necessary libraries
library(mvtnorm)
source('./scripts/eve-io.R')

# Parameter implementation
# Use this function to get the per gene parameter matrix
calculateParams <- function(gene.data, num.indivs)
{
  # Create the paramater matrix with the same number of rows as genes and one column for each parameter, and then name cols
  param.matrix <- matrix(nrow = nrow(gene.data), ncol = 4)
  colnames(param.matrix) <- c("theta", "sigma2", "alpha", "beta")
  
  # Create vectors in which to store the mean value of data per species per row and the variance within species
  species.mean <- vector(mode = "numeric", length = length(num.indivs))
  species.var <- vector(mode = "numeric", length = length(num.indivs))
  alpha <- .5
  
  # For each row within the paramater matrix, fill each column with parameter values
  for(row in 1:nrow(param.matrix))
  {
    # These start and end indices define where one species begins and ends for a gene
    start <- 1
    end <- num.indivs[1]
    
    for(i in 1:length(num.indivs))
    {
      # This increments the end point for a species for a gene
      if(i != 1)
      {
        end <- end + num.indivs[i]
      }
      # Calculate the mean of expression for each species and the variance of expression for each species
      species.mean[i] <- mean(gene.data[row, start:end])
      species.var[i] <- var(gene.data[row, start:end])
      # This increments the beginning point of a species for a gene
      start <- start + num.indivs[i]
    }
    # Fill the paramter matrix with theta, sigma sqaured, alpha, and beta values
    param.matrix[row, 1] <- mean(species.mean)
    param.matrix[row, 2] <- var(species.mean)
    param.matrix[row, 3] <- alpha
    param.matrix[row, 4] <- mean(species.var) / param.matrix[row, 2]
      
  }
  return(param.matrix)
}

# Expectation and variance from phylogeny
# Use this function to calculate expected species mean, and evolutionary
# variance for each node, given theta, alpha, and sigma^2 for each edge
calcExpVarOU <- function(tree, theta, alpha, sigma.sqaured, 
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
      sigma.sqaured / (2 * alpha) * (1 - exp(-2 * alpha * delta.T))
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
expandECovMatrix <- function(expected.mean, covar.matrix, sigma.sqaured, alpha, beta, num.indivs )
{
  index.expand <- rep(1:length(num.indivs), num.indivs)
  
  # Expand covariance matrix
  covar.matrix.exp <- covar.matrix[index.expand, index.expand]
  # Add within species variance
  diag(covar.matrix.exp) <- diag(covar.matrix.exp) + beta * sigma.sqaured / (2 * alpha)
  
  # Expand expected values
  expected.mean.exp <- expected.mean[index.expand]
  
  return(list(cov.matr = covar.matrix.exp, expected.mean = expected.mean.exp))
}

# Maximum likelihood estimations

logLikOU <- function(param.matrix.row, tree, gene.data.row, num.indivs)
{
  # Create vectors of paramters to pass into the function for calculating expression variance 
  theta <- param.matrix.row[1]
  sigma.squared <- param.matrix.row[2]
  alpha <- param.matrix.row[3]
  beta <- param.matrix.row[4]
  
  # Define the expected species mean for the root and the evolutionary variance for the root
  expected.species.mean.root <- theta
  evol.var.root <- sigma.squared / (2 * alpha)
  
  
  expression.var <- calcExpVarOU(tree, theta, alpha, sigma.squared, evol.var.root, expected.species.mean.root)
  covar.matrix <- calcCovMatOU(tree, alpha, expression.var$evol.variance)
  expanded.matrix <- expandECovMatrix(expression.var$expected.mean, covar.matrix, sigma.squared, alpha, beta, num.indivs)
  
  # Get log likelihood from the multivariate normal distribution density
  ll <- dmvnorm(x = gene.data.row, mean = expanded.matrix$expected.mean, sigma = expanded.matrix$cov.matr, log = TRUE )
  return(ll)
}

calculateLLPerGene <- function(param.matrix.row, tree, gene.data.row, num.indivs)
{
  ll <- -logLikOU(param.matrix.row, tree, gene.data.row, num.indivs)
  return(ll)
}

# Use this function to calculate the total likelihood 
calculateTotalLL <- function(ll.pergene)
{
  llTotal <- prod(ll.pergene)
  return(llTotal)
}

# Test for divergence and diversity between genes within species assuming an individual beta value for each gene
calculateLLIndivBeta <- function()
{
  #Initialize the tree
  tree <- initializePhylogeny()
  
  #Initialize the number of samples per species
  num.indivs <- getIndividuals()
  
  #Intialize the gene data
  gene.data <- getExprData(num.indivs)
  
  #Calculate the per gene parameter matrix based on the gene data
  param.matrix <- calculateParams(gene.data, num.indivs)
  
  # Create a vector to hold the values of the likelihood per gene and a list to hold the return values from the call to optim
  ll.pergene <- vector(mode = "numeric", length = nrow(param.matrix))
  max.params <- list()
  
  # For each gene, optimize the parameters and store the resulting likelihood in the likelihood vector
  for(row in 1:length(ll.pergene))
  {
    max.params <- optim(param.matrix[row, ], fn = calculateLLPerGene, gr = NULL, tree, gene.data[row, ], num.indivs,
                     method = "L-BFGS-B", lower = c(.001, .001))
    ll.pergene[row] <- as.numeric(max.params[2])
  }
  
  # Calculate the total likelihood as the product of all values within the likelihood vector
  ll.total <- calculateTotalLL(ll.pergene)
  
  # Save the results to a file that can be loaded into any future R environment for future use
  # To load this data simply call load("./results/llindivbetaresults.RData")
  save(ll.pergene, ll.total, file = "./results/llindivbetaresults.RData")
  
  return(ll.total)
}


