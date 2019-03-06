# This file contains functions to calculate gene likelihoods assuming a shared beta value between genes. 
# To run the calculation as a whole execute the function calculateLLsharedBeta. The total likelihood result
# will be returned from the function. The full results will be stored to an .RData file; to view these results 
# call load("./results/llsharedbetaresults.RData"). 
# Author: Rori Rohlfs, Lars Gronvold, John Mendoza
# Date 2/25/19

#Include all our necessary libraries
library(mvtnorm)
source('./scripts/eve-io.R')
source('./scripts/dvdt.R')

# This function calculates the initial parameters for the test where the beta paramater is assumed to be shared between genes 
calculateParamsSharedBeta <- function(gene.data, num.indivs)
{
  # Create the paramater matrix with the same number of rows as genes and one column for each parameter, and then name cols
  param.matrix <- matrix(nrow = nrow(gene.data), ncol = 4)
  colnames(param.matrix) <- c("theta", "sigma2", "alpha", "betaShared")
  
  # Create a vector in which to store the mean value of data per species per row
  species.mean <- vector(mode = "numeric", length = length(num.indivs))
  # Create a matrix to store the variance within species for each species and each gene
  species.var <- matrix(nrow = nrow(gene.data), ncol = length(num.indivs))
  # Create a vector to store the mean of the variance within species for each gene
  mean.species.var <- vector(mode = "numeric", length = nrow(species.var))
  
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
      species.var[row, i] <- var(gene.data[row, start:end])
      # This increments the beginning point of a species for a gene
      start <- start + num.indivs[i]
    }
    # Fill the paramter matrix with theta, sigma sqaured, and alpha values. Compute the mean variance for each gene and store it
    param.matrix[row, 1] <- mean(species.mean)
    param.matrix[row, 2] <- var(species.mean)
    param.matrix[row, 3] <- alpha
    mean.species.var[row] <- mean(species.var[row,])
    
  }
  
  # Compute and store the shared beta value as the mean of mean variance for each species for each gene divided by the mean of
  # sigma values
  param.matrix[,4] <- mean(mean.species.var) / mean(param.matrix[,2])
  return(param.matrix)
}

make.LogLikOU <- function(tree, gene.data.row, num.indivs, fixed = c(FALSE, FALSE, FALSE, FALSE))
{
  params <- fixed
  
  function(p)
  {
    # Use a boolean vector to determine which paramter to fix for the optimization
    params[!fixed] <- p[!fixed]
    
    # Create variables with paramters to pass into the function for calculating expression variance
    theta <- params[1]
    sigma.squared <- params[2]
    alpha <- params[3]
    beta <- params[4]

    # Define the expected species mean and evolutionary variance for the root
    expected.species.mean.root <- theta
    evol.var.root <- sigma.squared / (2 * alpha)
    
    expression.var <- calcExpVarOU(tree, theta, alpha, sigma.squared, evol.var.root, expected.species.mean.root)
    covar.matrix <- calcCovMatOU(tree, alpha, expression.var$evol.variance)
    expanded.matrix <- expandECovMatrix(expression.var$expected.mean, covar.matrix, sigma.squared, alpha, beta, num.indivs)
    
    # Get log likelihood from the multivariate normal distribution density
    ll <- dmvnorm(x = gene.data.row, mean = expanded.matrix$expected.mean, sigma = expanded.matrix$cov.matr, log = TRUE)
    return(-ll)
  }
}

# Test for divergence and diversity between genes within species assuming a shared beta value between genes
calculateLLsharedBeta <- function(tree, num.indivs, gene.data)
{
  # Calculate the inital per gene parameter matrix based on the gene data; note that in this version of the function, beta is not
  # calculated as we are under the assumption that beta is shared between genes
  param.matrix <- calculateParamsSharedBeta(gene.data, num.indivs)
  
  # Store the beta value in an independent variable and remove the beta column from the paramter matrix
  #beta <- param.matrix[1,4]
  #param.matrix <- param.matrix[, -4]
  
  # Create a vector to hold the values of the likelihood per gene and a list to hold the return values from the call to optim
  ll.pergene.sharedBeta <- vector(mode = "numeric", length = nrow(param.matrix))
  max.params <- list()
  
  # First, create an environment which contains the function to be optimized, and can hold all values other than beta fixed
  # during the optimization
  llOU <- make.LogLikOU(tree, gene.data[1, ], num.indivs, 
                        c(param.matrix[1, 1], param.matrix[1, 2], param.matrix[1, 3], FALSE))
  # Next, optimize the beta value for the first row of the initial paramters
  max.params <- optim(param.matrix[1, ], fn = llOU, gr = NULL, method = "L-BFGS-B", lower = c(.000000000000001))
  # Store the result in a numerical vector, and then set the beta value for all genes as this result
  vec <- as.numeric(max.params[[1]])
  param.matrix[, 4] <- vec[4]
  
  # For each gene, optimize the parameters keeping the shared beta fixed, and store the resulting likelihood in the likelihood vector
  for(row in 1:length(ll.pergene.sharedBeta))
  {
    # Now, create an enviroment with the function to be optimized, and hold the shared beta value fixed
    llOU <- make.LogLikOU(tree, gene.data[row, ], num.indivs, 
                          c(FALSE, FALSE, FALSE, param.matrix[row, 4]))
    # Store the results of the optimization in a list
    max.params <- optim(param.matrix[row, ], fn = llOU, gr = NULL, method = "L-BFGS-B", 
                        lower = c(.000000000000001, .000000000000001, .000000000000001))
    # Store the likelihood calculation based on the optimized parameter in a vector
    ll.pergene.sharedBeta[row] <- as.numeric(max.params[2])
  }
  
  # Calculate the total likelihood as the product of all values within the likelihood vector
  ll.total.sharedBeta <- calculateTotalLL(ll.pergene.sharedBeta)
  
  # Save the results to a file that can be loaded into any future R environment for future use
  # To load this data simply call load("./results/llsharedbetaresults.RData")
  save(ll.pergene.sharedBeta, ll.total.sharedBeta, file = "./results/llsharedbetaresults.RData")
  
  return(ll.total.sharedBeta)
}