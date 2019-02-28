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

# Test for divergence and diversity between genes within species assuming a shared beta value between genes
calculateLLsharedBeta <- function()
{
  # Initialize the tree
  tree <- initializePhylogeny()
  
  # Initialize the number of samples per species
  num.indivs <- getIndividuals()
  
  # Initialize the gene data
  gene.data <- getExprData(num.indivs)
  
  # Calculate the inital per gene parameter matrix based on the gene data; note that in this version of the function, beta is not
  # calculated as we are under the assumption that beta is shared between genes
  param.matrix <- calculateParamsSharedBeta(gene.data, num.indivs)
  
  # Store the beta value in an independent variable and remove the beta column from the paramter matrix
  #beta <- param.matrix[1,4]
  #param.matrix <- param.matrix[, -4]
  
  # Create a vector to hold the values of the likelihood per gene and a list to hold the return values from the call to optim
  ll.pergene.sharedBeta <- vector(mode = "numeric", length = nrow(param.matrix))
  max.params <- list()
  
  # For each gene, optimize the parameters and store the resulting likelihood in the likelihood vector
  for(row in 1:length(ll.pergene.sharedBeta))
  {
    max.params <- optim(param.matrix[row, ], fn = calculateLLPerGene, gr = NULL, tree, gene.data[row, ], num.indivs,
                        method = "L-BFGS-B", lower = c(.001, .001))
    ll.pergene.sharedBeta[row] <- as.numeric(max.params[2])
  }
  
  # Calculate the total likelihood as the product of all values within the likelihood vector
  ll.total.sharedBeta <- calculateTotalLL(ll.pergene.sharedBeta)
  
  # Save the results to a file that can be loaded into any future R environment for future use
  # To load this data simply call load("./results/llsharedbetaresults.RData")
  save(ll.pergene.sharedBeta, ll.total.sharedBeta, file = "./results/llsharedbetaresults.RData")
  
  return(ll.total.sharedBeta)
}