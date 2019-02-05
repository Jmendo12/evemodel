#Include all our necessary libraries
library(ape)

## Phylogengy Implementation
# Read the text in from a file name given by the user
fileName <- paste("./data/", sep = "", 
                  readline(prompt = "Enter the file in which the phylogeny is stored - this should be stored in a directory called data, and include the file extension: "))
# Create a tree from the file given by the user
tree <- read.tree(fileName)
# Reset the margins
par(mar = c(5, 4, 4, 2) + 0.1)
# Plot the tree so the edge, node numbers, and tip labels are shown
plot.phylo(tree, type = "p", show.tip.label = T, label.offset = .05,
           y.lim = c(0.5, 5.5), no.margin = F)
nodelabels(text = 1 : (Ntip(tree) + Nnode(tree)),
           node = 1 : (Ntip(tree) + Nnode(tree)), frame = "circle")
edgelabels(frame = "none", col = "red", adj = c(0.5, -0.1))
axisPhylo(1)

# Expectation and variance from phylogeny

# Use this function to calculate expected species mean (E), and evolutionary
# variance (V) for each node, given theta, alpha, and sigma^2 for each edge
calcExpVarOU <- function(tree, thetas, alphas, sigmas.sqaured, 
                         rootVar, rootE)
{
  # Declare vectors of expectation values and variances for all tips
  # and internal nodes
  E <- numeric( Ntip(tree) + Nnode(tree) )
  V <- numeric( Ntip(tree) + Nnode(tree) )
  
  # Get the order of the edges to traverse from root to tips
  edgeOrder <- rev(postorder(tree))
  
  # Get index of root
  root.index <- tree$edge[edgeOrder[1], 1]
  
  # Set root values
  E[root.index] = rootE
  V[root.index] = rootVar
  
  for(i in edgeOrder)
  {
    parent.index <- tree$edge[i, 1]
    child.index <- tree$edge[i, 2]
    delta.T <- tree$edge.length[i]
    E[child.index] <- E[parent.index] * exp(-alphas[i] * delta.T) +
      thetas[i] * (1 - exp(-alphas[i] * delta.T))
    V[child.index] <- V[parent.index] * exp(-2 * alphas[i] * delta.T) +
      sigmas.sqaured[i] / (2 * alphas[i]) * (1 - exp(-2 * alphas[i] * delta.T))
  }
  
  return(data.frame(E, V))
}

# Covaraince from phylogeny

# Use this function to calculate the covariance matrix given alpha for each edge
calcCovMatOU <- function(tree, alphas, V)
{
  # Copy the phylogeny and multiply edge lengths with alpha
  attenuationTree <- tree
  attenuationTree$edge.length <- attenuationTree$edge.length * alphas
  
  # calculate the attenuation matrix using the cophenetic distance function
  attenuation.Matrix <- cophenetic(attenuationTree)
  
  # get matrix of variances of the most recent common ancestors
  variance.MRCA <- apply(mrca(tree), 1:2, function(i) V[i])
  
  # Return the covariance matrix
  return (variance.MRCA * exp(-A))
}

# Use this function to add within-species variance. This will expand
# the covariance matrix and species mean vector to account for within-species variance
expandECovMatrix <- 