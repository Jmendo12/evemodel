# This file contains functions to calculate gene likelihoods assuming an individual beta value between genes.
# To run the calculation as a whole execute the function calculateLLIndivBeta. The total likelihood result
# will be returned from the function. The full results will be stored to an .RData file; to view these results
# call load("./results/llindivbetaresults.RData").
# Author: Rori Rohlfs, Lars Gronvold, John Mendoza
# Date 2/25/19


#' Likelihood for a gene assuming an individual beta.
#'
#' This function calculates the log likelihood for an individual gene assuming the gene has an individual beta value.
#'
#' Calculates the log likelihood for a gene with its own value for the beta parameter. The function calculates a log likelihood under
#' the OU model with four paramaters - theta, sigma squared, alpha, and beta. This function uses a model where each gene has
#' its own value for the beta parameter.
#'
#' @param param.matrix.row Numeric vector that contains the initial parameter estimates (theta, sigma2, alpha, beta)
#'     for a gene
#' @param tree Phylogenic tree created from a .newick file
#' @param gene.data.row Numeric vector that contains expression data of a gene
#' @param index.expand Numeric which specifies the index on which to expand a covariance matrix
#' @return A numeric value that represents the maximum likelihood of the gene.
#' @examples
#' optim(c(12.3, 6.1, .94, 3.6), fn = LLPerGeneIndivBeta,
#' tree, gene.data, index.expand)
#' optim(params, fn = LLPerGeneIndivBeta, tree, gene.data, index.expand)
LLPerGeneIndivBeta <- function(param.matrix.row, tree, gene.data.row, index.expand)
{
  ll <- -logLikOU( theta = param.matrix.row[1],
                   sigma2 = param.matrix.row[2],
                   alpha = param.matrix.row[3],
                   beta = param.matrix.row[4],
                   tree, gene.data.row, index.expand)
  return(ll)
}

#' Maximum likelihood estimation assuming an individual beta for each gene.
#'
#' Calculates the parameters to achieve the maximum likelihood for a gene assuming that each gene has it's own individual
#'     beta value.
#'
#' Calculates the maximum likelihood for a gene and its associated parameters through optimizing the original parameter
#'     estimates. This function assumes that each gene has its own individual beta value, and optimizes four parameters
#'     using optim, to find the parameters that give the greatest likelihood.
#'
#' @param tree a phylogenic tree created from a .newick file
#' @param gene.data a matrix which contains gene data
#' @param colSpecies a vector which specifies the species for columns of the gene data matrix. By default this parameter
#'     is set to the column names contained within the gene data matrix.
#' @return Returns a list which contains a list for each gene of their parameter values, the calculated maximum
#'     likelihood, the number of optim iterations, and a message on the exit status of the optim routine
#' @examples
#' test_res <- fitIndivBeta(tree, gene.data)
#' test_res <- fitIndivBeta(tree, gene.data, colSpecies = c("a", "b"))
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
      stats::optim(initial.param.matrix[row, ], fn = LLPerGeneIndivBeta, gr = NULL, tree, gene.data[row, ], index.expand,
            method = "L-BFGS-B", lower = c(-Inf, 1e-10, 1e-10, 1e-10), upper = c(Inf, Inf, alphaMax, Inf))
    }, error = function(e) {
      warning(paste(e$message, "at gene.data row", row), immediate. = T)
    })
  }) -> res

  return(res)
}

#' Likelihood for a gene assuming a shared beta value.
#'
#' This function calculates the log likelihood for genes assuming that each
#'     gene has the same beta value; a shared beta value between all genes.
#'
#'
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
    stats::optim( initial.param.matrix[row, 1:3], fn = LLPerGeneSharedBeta, method = "L-BFGS-B",
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
  sharedBetaFit <- stats::optimize(f = LLSharedBeta,interval=c(0.0001,100),
                            tree=tree, gene.data=gene.data)
  sharedBeta <- sharedBetaFit$minimum

  cat("fit with shared beta =",sharedBeta,"...\n")
  sharedBetaRes <- fitSharedBeta(sharedBeta, tree, gene.data, colSpecies)

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
