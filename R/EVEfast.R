prepEVEmodel <- function(tree, index.expand, thetaIdx, alphaIdx, sigma2Idx, betaIdx){
  
  # Force evaluation of arguments to avoid potential problems with lazy evaluation.
  force(index.expand); force(thetaIdx); force(alphaIdx); force(sigma2Idx); force(betaIdx)
  
  Nnodes <- Ntip(tree) + Nnode(tree)
  edgeOrder <- rev(postorder(tree))
  edge <- tree$edge
  edge.length <- tree$edge.length
  root.index <- edge[edgeOrder[1], 1]
  
  # Find the edge that defines the root regime!
  rootEdge <- edgeOrder[1]
  
  # Since regimes are defined per edge, there can be multiple regimes
  # associated with the root node. For now, one (random) of the edges from root define the root regime
  # TODO: Test if all root edges have same regime. Something like: 
  # if( length(unique(edgeRegimes[which(edge[,1]==root.index)])) > 1 ) stop("Regime change at root not supported!")

  # flattened matrix of precomputed mrca nodes
  mrcaV <- as.vector(mrca(tree))
  
  # edges to tip (used to calculate within-species variance)
  tipEdges <- match(1:Ntip(tree),edge[,2])
  
  # define EVEmodel function to return:
  function( par )
  {
    # vectorise parameters
    thetas <- par[thetaIdx]
    alphas <- par[alphaIdx]
    sigma2s <- par[sigma2Idx]
    beta <- par[betaIdx]
    
    ### 
    expected.mean <- numeric( Nnodes )
    evol.variance <- numeric( Nnodes )
    
    # Set root values
    expected.mean[root.index] = thetas[rootEdge]
    evol.variance[root.index] = sigma2s[rootEdge]/(2*alphas[rootEdge])
    
    # precalculate some expressions
    expAlphaT <- exp(-alphas * edge.length)
    expAlphaT2 <- exp(-2 * alphas * edge.length)
    tmp1 <- thetas * (1 - expAlphaT)
    tmp2 <- sigma2s / (2 * alphas) * (1 - expAlphaT2)

    for(i in edgeOrder)
    {
      parent.index <- edge[i, 1]
      child.index <- edge[i, 2]
      expected.mean[child.index] <- expected.mean[parent.index] * expAlphaT[i] + tmp1[i]
      evol.variance[child.index] <- evol.variance[parent.index] * expAlphaT2[i] + tmp2[i]
    }
    
    # Copy the phylogeny and multiply edge lengths with alphas
    attenuationTree <- tree
    attenuationTree$edge.length <- attenuationTree$edge.length * alphas
    
    # calculate the attenuation matrix using the cophenetic distance function
    attenuation.Matrix <- cophenetic(attenuationTree)
    
    # covariance matrix
    sigma <- evol.variance[mrcaV] * exp(-attenuation.Matrix)
    
    # Expand covariance matrix
    sigma <- sigma[index.expand, index.expand]
    # Add within-species variance
    diag(sigma) <- diag(sigma)*(1+beta) 

    # Return expanded sigma and mean (expected values)
    return(list(sigma=sigma, mean = expected.mean[index.expand]))
  }
}

dmvnorm_nocheck <- function(x, sigma, mean)
{
  if (is.vector(x)) 
    x <- matrix(x, ncol = length(x))
  p <- ncol(x)
  
  dec <- tryCatch(chol(sigma), error = function(e) e)
  if (inherits(dec, "error")) {
    x.is.mu <- colSums(t(x) != mean) == 0
    logretval <- rep.int(-Inf, nrow(x))
    logretval[x.is.mu] <- Inf
  }
  else {
    tmp <- backsolve(dec, t(x) - mean, transpose = TRUE)
    rss <- colSums(tmp^2)
    logretval <- -sum(log(diag(dec))) - 
      0.5 * p * log(2 * pi) - 0.5 * rss
  }

  logretval
}