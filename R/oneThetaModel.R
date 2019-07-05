#' Estimate initial parameter values for the one-theta model
#'
#' @param gene.data A matrix of expression values with samples in columns and genes in rows
#' @param colSpecies A character vector with same length as columns in gene.data, specifying the
#' species for the corresponding column.
#'
#' @return A matrix of initial parameters with the four parameters theta, sigma2, alpha and beta in columns and genes in rows
#' 
#' @importFrom stats var
initParamsOneTheta <- function(gene.data, colSpecies)
{
  colSpeciesIndices <- split(seq_along(colSpecies), f = colSpecies)
  species.mean <- sapply(colSpeciesIndices, function(i){ rowMeans(gene.data[,i]) })
  species.var <- sapply(colSpeciesIndices, function(i){ apply(gene.data[,i],1,var) })
  
  theta <- rowMeans(species.mean,na.rm = T)
  sigma2 <- apply(species.mean,1,var,na.rm = T)
  alpha <- .5
  beta <- rowMeans(species.var,na.rm = T) / sigma2
  
  return(cbind(theta,sigma2,alpha,beta))
}

#' Simulate expression data for the one-theta EVE model, i.e. same parameters on all edges
#'
#' @inheritParams simTwoTheta
#' @param theta Value of the theta parameter
#'
#' @return Matrix of simulated gene expression values with samples in columns and genes in rows
#' @export
#' @import mvtnorm
#' @importFrom stats setNames
simOneTheta <- function( n, tree, colSpecies, theta, sigma2, alpha, beta){
  Nedges <- Nedge(tree)
  
  index.expand <- match(colSpecies, tree$tip.label)
  
  mvdist <- EVEmodel(tree = tree, thetas = rep(theta,Nedges),alphas = rep(alpha,Nedges),sigma2s = rep(sigma2,Nedges),
                      beta = beta, index.expand = index.expand, rootE = theta, rootVar = beta * sigma2 / (2 * alpha))
  
  simData <- rmvnorm(n = n, mean = mvdist$mean, sigma = mvdist$sigma )
  colnames(simData) <- colSpecies
  return(simData)
}

#' Fit the maximum likelihood estimate of the parameters for a specific EVE model
#'
#' @param tree Phylogeny
#' @param gene.data A matrix of expression values with samples in columns and genes in rows
#' @param lowerBound A named numeric vector of the lower bound for the estimated parameters
#' @param upperBound A named numeric vector of the upper bound for the estimated parameters
#' @param logTransPars A character vector with the names of the parameters that will be log transformed
#' @param cores Number of parallel processes to run
#' @param fork Use forking for parallel execution
#' @param sharedBeta the shared beta parameter
#'
#' @return A list with:
#' \describe{
#'   \item{par}{A matrix with the parameter estimates}
#'   \item{ll}{The log likelihood for each estimate}
#'   \item{iterations}{Number of iterations of the optimisation routine before reaching the estimate}
#'   \item{convergence}{Error code from the \code{\link[stats]{optim}} function. 0 means convergence.}
#'   \item{message}{Error message from the \code{\link[stats]{optim}} function.}
#' }
#' @export
#' @import parallel
fitOneTheta <- function( tree, gene.data, colSpecies = colnames(gene.data), 
                         lowerBound = c(theta = -99, sigma2 = 0.0001, alpha = 0.001, beta = 0.001),
                         upperBound = c(theta =  99, sigma2 =   9999, alpha = 999  , beta = 99   ),
                         logTransPars = c("alpha","sigma2","beta"),
                         cores = 1, fork=F)
{
  #Calculate the per gene parameter matrix based on the gene data
  initPar <- initParamsOneTheta(gene.data, colSpecies)
  paramNames <- colnames(initPar)
  
  
  # log transform initial parameters and bounds
  doTransPar <- paramNames %in% logTransPars # which params to log transform
  initPar[,doTransPar] <- log(initPar[,doTransPar])
  lowerBound[logTransPars] <- log(lowerBound[logTransPars])
  upperBound[logTransPars] <- log(upperBound[logTransPars])
  
  # match the column species with the phylogeny tip labels
  index.expand <- match(colSpecies, tree$tip.label)
  
  localEVEmodel <- prepEVEmodel(tree = tree,index.expand = index.expand,
                                thetaIdx = rep(match("theta",paramNames),Nedge(tree)),
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
