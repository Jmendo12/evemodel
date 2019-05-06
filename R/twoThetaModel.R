#' Estimate initial parameter values for two-theta model without using phylogeny
#'
#' @param gene.data A matrix of expression values with samples in columns and genes in rows
#' @param colSpecies A character vector with same length as columns in gene.data, specifying the
#' species for the corresponding column.
#' @param shiftSpecies character vector with species that are assigned to theta2
#'
#' @return A matrix of initial parameters with the parameters in columns and genes in rows
initParamsTwoTheta <- function(gene.data, colSpecies, shiftSpecies)
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

simTwoTheta <- function( n, tree, colSpecies, isTheta2edge, theta1, theta2, sigma2, alpha, beta){
  Nedges <- Nedge(tree)
  
  index.expand <- match(colSpecies, tree$tip.label)
  
  # TODO: check if theta2 is on root edge. Now it assumes that root is theta1
  mvdist <- EVEmodel(tree = tree, thetas = ifelse(isTheta2edge, theta2,theta1),
                     alphas = rep(alpha,Nedges),sigma2s = rep(sigma2,Nedges),
                     beta = beta, index.expand = index.expand, rootE = theta1, rootVar = beta * sigma2 / (2 * alpha))
  
  simData <- rmvnorm(n = 100, mean = mvdist$mean, sigma = mvdist$sigma )
  colnames(simData) <- colSpecies
  return(simData)
}

fitTwoTheta <- function( tree, gene.data, isTheta2edge, colSpecies = colnames(gene.data), 
                         lowerBound = c(theta1 = -99, theta2 = -99, sigma2 = 0.0001, alpha = 0.001, beta = 0.001),
                         upperBound = c(theta1 =  99, theta2 =  99, sigma2 =   9999, alpha = 999  , beta = 99   ),
                         logTransPars = c("alpha","sigma2","beta"),
                         cores = 1, fork=F)
{
  # get shift species from shift edges
  shiftSpecies = tree$tip.label[tree$edge[isTheta2edge & tree$edge[,2] <= Ntip(tree),2]]
  
  initPar <- initParamsTwoTheta(gene.data, colSpecies, shiftSpecies = shiftSpecies)
  paramNames <- colnames(initPar)
  
  
  # log transform initial parameters and bounds
  doTransPar <- paramNames %in% logTransPars # which params to log transform
  initPar[,doTransPar] <- log(initPar[,doTransPar])
  lowerBound[logTransPars] <- log(lowerBound[logTransPars])
  upperBound[logTransPars] <- log(upperBound[logTransPars])
  
  # match the column species with the phylogeny tip labels
  index.expand <- match(colSpecies, tree$tip.label)
  
  localEVEmodel <- prepEVEmodel(tree = tree,index.expand = index.expand,
                                thetaIdx = ifelse(isTheta2edge, match("theta2",paramNames),
                                                  match("theta1",paramNames)),
                                alphaIdx = rep(match("alpha",paramNames),Nedge(tree)),
                                sigma2Idx = rep(match("sigma2",paramNames),Nedge(tree)),
                                betaIdx = match("beta",paramNames))
  
  fitOneGene <- function(row){
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
  }
  
  myFunc <- function(row) {localEVEmodel(par = initPar[row,])}
  
  
  if(cores==1){
    res <- lapply(X = 1:nrow(gene.data), FUN = fitOneGene)
  } else if(fork){
    res <- mclapply(mc.cores = cores,X = 1:nrow(gene.data), FUN = fitOneGene)
  }else{
    cl <- makeCluster(cores)
    # Export local environment to worker processes
    clusterExport(cl, varlist = c(ls(envir = environment()),"dmvnorm_nocheck"),envir = environment())
    # clusterExport(cl, varlist = c("initPar","gene.data","lowerBound","upperBound","doTransPar","localEVEmodel","dmvnorm_nocheck"),envir = environment())
    clusterEvalQ(cl, expr = library(ape))
    res <- parLapply(cl = cl,X = 1:nrow(gene.data), fun = fitOneGene)
    stopCluster(cl)
  }
  
  # Simplify the results
  list( par = t(sapply(res,function(x) x$par)), 
        ll = -sapply(res,function(x) x$value),
        iterations = setNames(sapply(res,function(x) x$counts[1]),NULL),
        convergence = sapply(res,function(x) x$convergence),
        message = sapply(res,function(x) x$message))
}