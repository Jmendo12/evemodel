#Include all our necessary libraries
library(ape)

## Phylogengy Implementation
# Initiliaze the phylogeny data and plot the phylogeny
initializePhylogeny <- function()
{
  # Read in the filename the phylogeny is stored in
  fileName <- paste("./data/", sep = "", readline(prompt = "Enter the file in which the phylogeny data is stored, and include the file extension: "))
  # Create a tree from the file given by the user
  tree <- read.tree(fileName)
  # Reset the margins
  par(mar = c(5, 4, 4, 2) + 0.1)
  # Plot the tree so the edge, node numbers, and tip labels are shown
  plot.phylo(tree, type = "p", show.tip.label = T, label.offset = .05, y.lim = c(.5, 5.5), no.margin = F)
  nodelabels(text = 1 : (Ntip(tree) + Nnode(tree)), node = 1 : (Ntip(tree) + Nnode(tree)), frame = "circle")
  edgelabels(frame = "none", col = "red", adj = c(.5, -.1))
  axisPhylo(1)
  
  return(tree)
}

# Use this function to get the data on the number of samples per species
getIndividuals <- function()
{
  fileName <- paste("./data/", sep = "", readline(prompt = "Enter the file in which the number of samples per species data is stored, and include the file extension: "))
  num.indivs <- scan(fileName, sep = " ")
  return(num.indivs)
}

# Use this function to get the per gene expression data
getExprData <- function(num.indivs)
{
  fileName <- paste("./data/", sep = "", readline(prompt = "Enter the file in which the gene expression data is stored, and include the file extension: "))
  gene.data <- as.matrix(read.table("data/sampleExpr.dat",skip = 1,header = F,row.names = 1))
  
  colnames(gene.data) <- rep(LETTERS[1:5], num.indivs)
  
  return(gene.data)
}

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
  count <- 1
  
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
    param.matrix[row, 3] <- .5
    param.matrix[row, 4] <- mean(species.var) / param.matrix[row, 2]
      
  }
  return(param.matrix)
}

# Expectation and variance from phylogeny
# Use this function to calculate expected species mean, and evolutionary
# variance for each node, given theta, alpha, and sigma^2 for each edge
calcExpVarOU <- function(tree, thetas, alphas, sigmas.sqaured, 
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
      sigmas.sqaured[i] / (2 * alphas[i]) * (1 - exp(-2 * alphas[i] * delta.T))
  }
  
  return(data.frame(expected.mean, evol.variance))
}

# Covaraince from phylogeny

# Use this function to calculate the covariance matrix given alpha for each edge
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
  return (variance.MRCA * exp(-A))
}

# Use this function to add within-species variance. This will expand
# the covariance matrix and species mean vector to account for within-species variance
expandECovMatrix <- function(expected.mean, covar.matrix, sigma.sqaured, alpha, beta, num.indiv )
{
  index.expand <- rep(1:length(num.indiv), num.indiv)
  
  # Expand covariance matrix
  covar.matrix <- covar.matrix[index.expand, index.expand]
  # Add within species variance
  diag(covar.matrix) <- diag(covar.matrix) + beta * sigma.sqaured / (2 * alpha)
  
  # Expand expected values
  expected.mean <- expected.mean[index.expand]
  
  return(list(covar.matrix, expected.mean))
}

divergence.diversity.test <- function()
{
  #Initialize the tree
  tree <- initializePhylogeny()
  
  #Initialize the number of samples per species
  num.indivs <- getIndividuals()
  
  #Intialize the per gene data
  gene.data <- getExprData(num.indivs)
  
  #Calculate the per gene parameter matrix based on the gene data
  paramater.matrix <- calculate.params(gene.data)
}