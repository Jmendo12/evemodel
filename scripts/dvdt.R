# This file contains functions to calculate gene likelihoods assuming an individual beta value between genes. 
# To run the calculation as a whole execute the function calculateLLIndivBeta. The total likelihood result
# will be returned from the function. The full results will be stored to an .RData file; to view these results 
# call load("./results/llindivbetaresults.RData"). 
# Author: Rori Rohlfs, Lars Gronvold, John Mendoza
# Date 2/25/19

#Include all our necessary libraries
source('scripts/EVEcore.R')

LLPerGeneIndivBeta <- function(param.matrix.row, tree, gene.data.row, index.expand)
{
  ll <- -logLikOU(param.matrix.row, tree, gene.data.row, index.expand)
  return(ll)
}


# Maximum liklihood estimation of parameters with individual betas per gene
fitIndivBeta <- function(tree, gene.data, colSpecies = colnames(gene.data))
{
  #Calculate the per gene parameter matrix based on the gene data
  initial.param.matrix <- initialParams(gene.data, colSpecies)
  
  # match the column species with the phylogeny tip labels
  index.expand <- match(colSpecies, tree$tip.label)
  # Calculate max alpha based on the proportion of variance from covariance
  alphaMax <- -log(.01) / min(tree$edge.length[tree$edge[,2] <= Ntip(tree)])
  
  # For each gene, optimize the parameters and store the resulting likelihood in the likelihood vector
  lapply(1:nrow(gene.data), function(row){
    # Error handling to catch infinte optim or function values that arise when data with NaN paramters is optimized
    res <- tryCatch({
      optim(initial.param.matrix[row, ], fn = LLPerGeneIndivBeta, gr = NULL, tree, gene.data[row, ], index.expand,
            method = "L-BFGS-B", lower = c(-Inf, 1e-10, 1e-10, 1e-10), upper = c(Inf, Inf, alphaMax, Inf))
    }, error = function(e) {
      warning(paste(e$message, "at gene.data row", row), immediate. = T)
    })
  }) -> res

  return(res)
}


