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
  ll <- -logLikOU( theta = param.matrix.row[1],
                   sigma2 = param.matrix.row[2],
                   alpha = param.matrix.row[3],
                   beta = param.matrix.row[4],
                   tree, gene.data.row, index.expand)
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

LLPerGeneSharedBeta <- function(param.matrix.row, betaShared, tree, gene.data.row, index.expand)
{
  ll <- logLikOU( theta = param.matrix.row[1],
                  sigma2 = param.matrix.row[2],
                  alpha = param.matrix.row[3],
                  beta = betaShared,
                  tree, gene.data.row, index.expand)
  return(-ll)
}

fitSharedBeta <- function(betaShared, tree, gene.data, colSpecies = colnames(gene.data))
{
  #Calculate the per gene parameter matrix based on the gene data
  initial.param.matrix <- initialParams(gene.data, colSpecies)
  
  # match the column species with the phylogeny tip labels
  index.expand <- match(colSpecies, tree$tip.label)
  # Calculate max alpha based on the proportion of variance from covariance
  alphaMax <- -log(.01) / min(tree$edge.length[tree$edge[,2] <= Ntip(tree)])
  
  
  # For each gene, optimize the parameters and store the resulting likelihood in the likelihood vector
  lapply(1:nrow(gene.data), function(row){
    optim( initial.param.matrix[row, 1:3], fn = LLPerGeneSharedBeta, method = "L-BFGS-B",
           lower = c(-Inf, 1e-10, 1e-10, 1e-10), upper = c(Inf, Inf, alphaMax, Inf),
           betaShared= betaShared,
           tree = tree, gene.data.row = gene.data[row, ], index.expand=index.expand)
  }) -> res
  
  return(res)
}

LLSharedBeta <- function(betaShared, ...)
{
  cat("LLSharedBeta: beta =",betaShared) 

  resSharedBeta <- fitSharedBeta(betaShared, ...)
  
  sumLL <- sum(sapply(resSharedBeta, function(res) res$value ))
  
  nNotConverged <- sum(sapply(resSharedBeta, function(res) res$convergence )!=0)
  
  if( nNotConverged>0 ){
    cat("  ",nNotConverged,"gene(s) did not converge!") 
  }
  
  cat("  LL =",sumLL,"\n") 
  
  # return sum of -LL for all genes
  return(sumLL)
}


betaSharedTest <- function(tree, gene.data, colSpecies = colnames(gene.data)){
  cat("fit with individual betas...\n")
  indivBetaRes <- fitIndivBeta(tree,gene.data,colSpecies)
  
  cat("Estimate shared beta...\n")
  sharedBetaFit <- optimize(f = LLSharedBeta,interval=c(0.0001,100),
                            tree=tree, gene.data=gene.data)
  sharedBeta <- sharedBetaFit$minimum
  
  cat("fit with shared beta =",sharedBeta,"...\n")
  sharedBetaRes <- fitSharedBeta(sharedBeta, tree, gene.data)
  
  # calculate likelihood ratio test statistic
  LRT <- mapply(function(indivBetaRow, sharedBetaRow)
    {
      (2 * -(indivBetaRow$value)) - (2 * -(sharedBetaRow$value))
    }, indivBetaRes, sharedBetaRes)

  
  return( list(indivBetaRes = indivBetaRes, 
               sharedBeta = sharedBeta, 
               sharedBetaRes = sharedBetaRes, 
               LRT = LRT) )
}
