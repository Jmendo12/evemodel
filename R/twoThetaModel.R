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
  species.mean <- sapply(colSpeciesIndices, function(i){ rowMeans(gene.data[,i, drop=F]) })
  species.var <- sapply(colSpeciesIndices, function(i){ apply(gene.data[,i, drop=F],1,var) })
  
  
  nonShiftSpecies <- setdiff(colnames(species.mean),shiftSpecies)
  
  theta1 <- rowMeans(species.mean[ ,nonShiftSpecies, drop=F],na.rm=T)
  theta2 <- rowMeans(species.mean[ ,shiftSpecies, drop=F],na.rm=T)
  sigma2 <- apply(species.mean,1,var,na.rm=T)
  alpha <- .5
  beta <- rowMeans(species.var,na.rm=T) / sigma2
  
  return(cbind(theta1,theta2,sigma2,alpha,beta))
}

#' Simulate expression data for the two-theta EVE model, i.e. two different thetas assigned to specific edges
#'
#' @param n Number of "genes" to simulate
#' @param tree Phylogeny
#' @param colSpecies A character vector with same length as columns in returned expression matrix, 
#' specifying the species, i.e. tip labels in the phylogeny, for the corresponding column.
#' @param isTheta2edge Logical vector with same length as number of edges in tree specifying whether
#' the corresponding edge theta parameter should be theta2 (TRUE) or theta1 (FALSE)
#' @param theta1 Value of the theta1 parameter
#' @param theta2 Value of the theta2 parameter
#' @param sigma2 Value of the sigma2 parameter
#' @param alpha Value of the alpha parameter
#' @param beta Value of the beta parameter
#'
#' @return Matrix of simulated gene expression values with samples in columns and genes in rows
#' @export
simTwoTheta <- function( n, tree, colSpecies, isTheta2edge, theta1, theta2, sigma2, alpha, beta){
  Nedges <- Nedge(tree)
  
  index.expand <- match(colSpecies, tree$tip.label)
  
  # TODO: check if theta2 is on root edge. Now it assumes that root is theta1
  mvdist <- EVEmodel(tree = tree, thetas = ifelse(isTheta2edge, theta2,theta1),
                     alphas = rep(alpha,Nedges),sigma2s = rep(sigma2,Nedges),
                     beta = beta, index.expand = index.expand, rootE = theta1, rootVar = beta * sigma2 / (2 * alpha))
  
  simData <- rmvnorm(n = n, mean = mvdist$mean, sigma = mvdist$sigma )
  colnames(simData) <- colSpecies
  return(simData)
}


#' @describeIn fitOneTheta Fit model with two different thetas assigned to specific edges
#' @param isTheta2edge Logical vector with same length as number of edges in tree specifying whether
#' the corresponding edge theta parameter should be theta2 (TRUE) or theta1 (FALSE)
#' @param colSpecies A character vector with same length as columns in returned expression matrix, 
#' specifying the species, i.e. tip labels in the phylogeny, for the corresponding column.
#' @export
fitTwoTheta <- function( tree, gene.data, isTheta2edge, colSpecies = colnames(gene.data), 
                         extra.var = NULL,
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
                   extra.var.row = if(is.null(extra.var)) NULL else extra.var[row, ],
                   fn = function(par, gene.data.row, extra.var.row){
                     # reverse log transform parameters
                     par[doTransPar] <- exp(par[doTransPar])
                     
                     mvnormParams <- localEVEmodel(par)
                     
                     # Add extra variance (if given)
                     if( !is.null(extra.var) )
                       diag(mvnormParams$sigma) <- diag(mvnormParams$sigma) +  extra.var.row
                     
                     # ignore species with NA in the expression matrix
                     notNA <- !is.na(gene.data.row)
                     return(-dmvnorm_nocheck(gene.data.row[notNA], sigma = mvnormParams$sigma[notNA,notNA], 
                                             mean=mvnormParams$mean[notNA]))
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