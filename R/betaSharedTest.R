# This file contains functions to calculate gene likelihoods assuming an individual beta value between genes.
# To run the calculation as a whole execute the function calculateLLIndivBeta. The total likelihood result
# will be returned from the function. The full results will be stored to an .RData file; to view these results
# call load("./results/llindivbetaresults.RData").
# Author: Rori Rohlfs, Lars Gronvold, John Mendoza
# Date 2/25/19



#' @describeIn fitOneTheta Fit model with a given shared beta for all genes
#' @export
fitSharedBeta <- function( sharedBeta, tree, gene.data, colSpecies = colnames(gene.data), 
                           extra.var = NULL,
                         lowerBound = c(theta = -99, sigma2 = 0.0001, alpha = 0.001),
                         upperBound = c(theta =  99, sigma2 =   9999, alpha = 999 ),
                         logTransPars = c("alpha","sigma2","beta"),
                         cores = 1, fork=F)
{
  #Calculate the per gene parameter matrix based on the gene data
  initPar <- initParamsOneTheta(gene.data, colSpecies)[,c("theta","sigma2","alpha")]
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
                                betaIdx = 4)
  
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
                     
                     mvnormParams <- localEVEmodel(c(par, sharedBeta))
                     
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





#' Beta shared test
#'
#' @param tree Species phylogeny (phylo object)
#' @param gene.data A matrix of expression values with samples in columns and genes in rows
#' @param colSpecies A character vector with same length as columns in gene.data, specifying the
#' species for the corresponding column.
#' @param ... Parameters passed to \code{\link{fitOneTheta}} and \code{\link{fitSharedBeta}}
#'
#' @return List with:
#' \itemize{
#'   \item indivBetaRes: results from \code{\link{fitOneTheta}}
#'   \item sharedBetaRes: results from \code{\link{fitSharedBeta}}
#'   \item sharedBeta: the estimated shared beta parameter
#'   \item LRT: log likelihood ratio test statistic between the individual and shared beta model
#' }
#' @export
betaSharedTest <- function(tree, gene.data, colSpecies = colnames(gene.data), ...){
  cat("fit with individual betas...\n")
  indivBetaRes <- fitOneTheta(tree,gene.data,colSpecies, ...)

  LLSharedBeta <- function(betaShared, ...)
  {
    cat("LLSharedBeta: beta =",betaShared)
    
    resSharedBeta <- fitSharedBeta(betaShared, ...)
    
    sumLL <- sum(resSharedBeta$ll)
    
    nNotConverged <- sum(resSharedBeta$convergence!=0)
    
    if( nNotConverged>0 ){
      cat("  ",nNotConverged,"gene(s) did not converge!")
    }
    
    cat("  LL =",sumLL,"\n")
    
    # return -sum of LL for all genes
    return(-sumLL)
  }
  
  cat("Estimate shared beta...\n")
  sharedBetaFit <- stats::optimize(f = LLSharedBeta,interval=c(0.0001,100),
                            tree=tree, gene.data=gene.data, colSpecies=colSpecies, ...)
  sharedBeta <- sharedBetaFit$minimum

  cat("fit with shared beta =",sharedBeta,"...\n")
  sharedBetaRes <- fitSharedBeta(sharedBeta, tree, gene.data, colSpecies, ...)

  # calculate likelihood ratio test statistic
  LRT <- 2 * (indivBetaRes$ll - sharedBetaRes$ll)


  return( list(indivBetaRes = indivBetaRes,
               sharedBetaRes = sharedBetaRes,
               sharedBeta = sharedBeta,
               LRT = LRT) )
}
