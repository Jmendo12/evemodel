prepEVEmodel <- function(tree, index.expand, thetaIdx, alphaIdx, sigma2Idx, betaIdx){
  
  # Force evaluation of arguments. This helps detect syntax errors in the arguments earlier
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
    
    # for(i in edgeOrder)
    # {
    #   parent.index <- edge[i, 1]
    #   child.index <- edge[i, 2]
    #   delta.T <- edge.length[i]
    #   expected.mean[child.index] <- expected.mean[parent.index] * exp(-alphas[i] * delta.T) +
    #     thetas[i] * (1 - exp(-alphas[i] * delta.T))
    #   evol.variance[child.index] <- evol.variance[parent.index] * exp(-2 * alphas[i] * delta.T) +
    #     sigma2s[i] / (2 * alphas[i]) * (1 - exp(-2 * alphas[i] * delta.T))
    # }
    
    
    # precalculate some expressions
    expAlphaT <- exp(-alphas * edge.length)
    # expAlphaT2 <- expAlphaT^2  # This introduces some numerical differences, but is slightly faster
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
    # diag(sigma) <- diag(sigma) + beta * sigma2s[tipEdges] / (2 * alphas[tipEdges]) # tipEdges should have been expanded
    diag(sigma) <- diag(sigma)*(1+beta) # is this correct? It gives slight numerical differences that can result in diferrent paramet estimates

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

fitIndivBeta_fast <- function(tree, gene.data, colSpecies = colnames(gene.data))
{
  #Calculate the per gene parameter matrix based on the gene data
  initPar <- initialParams(gene.data, colSpecies)
  params <- colnames(initPar)
  
  # match the column species with the phylogeny tip labels
  index.expand <- match(colSpecies, tree$tip.label)
  
  localEVEmodel <- prepEVEmodel(tree = tree,index.expand = index.expand,
                                thetaIdx = rep(match("theta",params),Nedge(tree)),
                                alphaIdx = rep(match("alpha",params),Nedge(tree)),
                                sigma2Idx = rep(match("sigma2",params),Nedge(tree)),
                                betaIdx = match("beta",params))
  
  # Calculate max alpha based on the proportion of variance from covariance
  alphaMax <- -log(.01) / min(tree$edge.length[tree$edge[,2] <= Ntip(tree)])
  
  # For each gene, optimize the parameters and store the resulting likelihood in the likelihood vector
  lapply(1:nrow(gene.data), function(row){
    # Error handling to catch infinte optim or function values that arise when data with NaN paramters is optimized
    res <- tryCatch({
      stats::optim(par = initPar[row, ], 
                   method = "L-BFGS-B", lower = c(-Inf, 1e-10, 1e-10, 1e-10), 
                   upper = c(Inf, Inf, alphaMax, Inf),
                   gene.data.row = gene.data[row, ],
                   fn = function(par, gene.data.row){
                     mvnormParams <- localEVEmodel(par)
                     return(-dmvnorm_nocheck(gene.data.row, sigma = mvnormParams$sigma, mean=mvnormParams$mean))
                   })
    }, error = function(e) {
      warning(paste(e$message, "at gene.data row", row), immediate. = T)
    })
  }) -> res
  
  return(res)
}


# Include log transform and bounds
fitIndivBeta_fast_log <- function(tree, gene.data, colSpecies = colnames(gene.data), 
                                  lowerBound = c(theta = -99, sigma2 = 0.0001, alpha = 0.001, beta = 0.001),
                                  upperBound = c(theta =  99, sigma2 =   9999, alpha = 999  , beta = 99   ))
{
  #Calculate the per gene parameter matrix based on the gene data
  initPar <- initialParams(gene.data, colSpecies)
  paramNames <- colnames(initPar)
  
  
  # log transform initial parameters and bounds
  doTransPar <- paramNames %in% c("alpha","sigma2","beta") # which params to log transform
  initPar[,doTransPar] <- log(initPar[,doTransPar])
  lowerBound[c("alpha","sigma2","beta")] <- log(lowerBound[c("alpha","sigma2","beta")])
  upperBound[c("alpha","sigma2","beta")] <- log(upperBound[c("alpha","sigma2","beta")])
  
  # match the column species with the phylogeny tip labels
  index.expand <- match(colSpecies, tree$tip.label)
  
  localEVEmodel <- prepEVEmodel(tree = tree,index.expand = index.expand,
                               thetaIdx = rep(match("theta",paramNames),Nedge(tree)),
                               alphaIdx = rep(match("alpha",paramNames),Nedge(tree)),
                               sigma2Idx = rep(match("sigma2",paramNames),Nedge(tree)),
                               betaIdx = match("beta",paramNames))
  
  # For each gene, optimize the parameters and store the resulting likelihood in the likelihood vector
  lapply(1:nrow(gene.data), function(row){
    # Error handling to catch infinte optim or function values that arise when data with NaN paramters is optimized
    res <- tryCatch({
      stats::optim(par = initPar[row, ], method = "L-BFGS-B", 
                   lower = lowerBound, upper = upperBound,
                   gene.data.row = gene.data[row, ],
                   fn = function(par, gene.data.row){
                     # reverse log transform parameters
                     par[doTransPar] <- exp(par[doTransPar])
                     
                     mvnormParams <- localEVEmodel(par)
                     return(-dmvnorm_nocheck(gene.data.row, sigma = mvnormParams$sigma, mean=mvnormParams$mean))
                   })
    }, error = function(e) {
      warning(paste(e$message, "at gene.data row", row), immediate. = T)
    })
    # reverse log transform estimated parameters
    res$par[doTransPar] <- exp(res$par[doTransPar])
    return(res)
  }) -> res
  
  # TODO: convert to data.frame?
  # or maybe list of matrices?
  
  
  return(res)
}

