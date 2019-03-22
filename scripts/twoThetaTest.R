getEdgesFromMRCA <- function(tree, tips){
  mrcaNode <- getMRCA(tree,tips)

  getDescendantEdgesRecursive <- function(nodes){
    descendantEdges <- which(tree$edge[ ,1] %in% nodes)
    if( length(descendantEdges) > 0){
      descendantEdges <- c(descendantEdges, Recall(nodes = tree$edge[descendantEdges,2]))
    }
    return(descendantEdges)
  }
  
  return( getDescendantEdgesRecursive(mrcaNode) )
}

#' Estimate initial parameter values for two-theta model without using phylogeny
#'
#' @param gene.data A matrix of expression values with samples in columns and genes in rows
#' @param colSpecies A character vector with same length as columns in gene.data, specifying the 
#' species for the corresponding column.
#' @param shiftSpecies character vector with species that are assigned to theta2
#'
#' @return A matrix of initial parameters with the parameters in columns and genes in rows
initialParamsTwoTheta <- function(gene.data, colSpecies, shiftSpecies)
{
  colSpeciesIndices <- split(seq_along(colSpecies), f = colSpecies)
  species.mean <- sapply(colSpeciesIndices, function(i){ rowMeans(gene.data[,i]) })
  species.var <- sapply(colSpeciesIndices, function(i){ apply(gene.data[,i],1,var) })
  
  
  nonShiftSpecies <- setdiff(colnames(species.mean),shiftSpecies)
  
  theta1 <- rowMeans(species.mean[ ,nonShiftSpecies])
  theta2 <- rowMeans(species.mean[ ,shiftSpecies])
  sigma2 <- apply(species.mean,1,var)
  alpha <- .5
  beta <- rowMeans(species.var) / sigma2
  
  return(cbind(theta1,theta2,sigma2,alpha,beta))
}

fitTwoThetas <- function(tree, gene.data, shiftSpecies, colSpecies){
  # Calculate initial estimate of parameters
  initial.param.matrix <- initialParamsTwoTheta(gene.data, colSpecies, shiftSpecies)
  
  isThetaShiftEdge <- 1:Nedge(tree) %in% getEdgesFromMRCA(tree, tips = shiftSpecies)
  
  # match the column species with the phylogeny tip labels
  index.expand <- match(colSpecies, tree$tip.label)
  
  # Calculate max alpha based on the proportion of variance from covariance
  alphaMax <- -log(.01) / min(tree$edge.length[tree$edge[,2] <= Ntip(tree)])

  # function to optimize
  LLPerGeneTwoTheta <- function(params, gene.data.row)
  {
    ll <- logLikTwoTheta( theta1 = params[1], theta2 = params[2], sigma2 = params[3],
                          alpha = params[4], beta = params[5],
                          tree = tree, isThetaShiftEdge = isThetaShiftEdge, 
                          gene.data.row = gene.data.row, index.expand = index.expand)
    return(-ll)
  }
    
  # For each gene, optimize the parameters
  lapply(1:nrow(gene.data), function(row){
    optim( initial.param.matrix[row, ], fn = LLPerGeneTwoTheta, method = "L-BFGS-B",
           lower = c(-Inf, -Inf, 1e-10, 1e-10, 1e-10), upper = c(Inf, Inf, Inf, alphaMax, Inf),
           gene.data.row = gene.data[row, ])
  }) -> res
  
  return(res)
}


# test for shift in theta at branch specified by the extant species "shiftSpecies"
twoThetaTest <- function(tree, gene.data, shiftSpecies, colSpecies = colnames(gene.data)){
  
  cat("fit with single theta...\n") # same as individual beta
  oneThetaRes <- fitIndivBeta(tree,gene.data,colSpecies)
  
  cat("Fit with two thetas\n")
  twoThetaRes <- fitTwoThetas(tree, gene.data, shiftSpecies, colSpecies)
  
  # calculate likelihood ratio test statistic
  oneThetaLL <- -sapply(oneThetaRes, function(res) res$value )
  twoThetaLL <- -sapply(twoThetaRes, function(res) res$value )
  LRT <- 2 * (twoThetaLL - oneThetaLL)
  
  return( list(oneThetaRes = oneThetaRes, 
               twoThetaRes = twoThetaRes, 
               LRT = LRT) )
}
