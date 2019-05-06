# This file contains functions to calculate gene likelihoods assuming an individual beta value between genes.
# To run the calculation as a whole execute the function calculateLLIndivBeta. The total likelihood result
# will be returned from the function. The full results will be stored to an .RData file; to view these results
# call load("./results/llindivbetaresults.RData").
# Author: Rori Rohlfs, Lars Gronvold, John Mendoza
# Date 2/25/19


fitSharedBeta <- function( sharedBeta, tree, gene.data, colSpecies = colnames(gene.data), 
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
                   fn = function(par, gene.data.row){
                     # reverse log transform parameters
                     par[doTransPar] <- exp(par[doTransPar])
                     
                     mvnormParams <- localEVEmodel(c(par, sharedBeta))
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

# fitSharedBeta <- function(betaShared, tree, gene.data, colSpecies = colnames(gene.data))
# {
#   #Calculate the per gene parameter matrix based on the gene data
#   initial.param.matrix <- initialParams(gene.data, colSpecies)
# 
#   # match the column species with the phylogeny tip labels
#   index.expand <- match(colSpecies, tree$tip.label)
#   # Calculate max alpha based on the proportion of variance from covariance
#   alphaMax <- -log(.01) / min(tree$edge.length[tree$edge[,2] <= Ntip(tree)])
# 
# 
#   # For each gene, optimize the parameters and store the resulting likelihood in the likelihood vector
#   lapply(1:nrow(gene.data), function(row){
#     stats::optim( initial.param.matrix[row, 1:3], fn = LLPerGeneSharedBeta, method = "L-BFGS-B",
#            lower = c(-Inf, 1e-10, 1e-10, 1e-10), upper = c(Inf, Inf, alphaMax, Inf),
#            betaShared= betaShared,
#            tree = tree, gene.data.row = gene.data[row, ], index.expand=index.expand)
#   }) -> res
# 
#   return(res)
# }

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


betaSharedTest <- function(tree, gene.data, colSpecies = colnames(gene.data)){
  cat("fit with individual betas...\n")
  indivBetaRes <- fitOneTheta(tree,gene.data,colSpecies)

  cat("Estimate shared beta...\n")
  sharedBetaFit <- stats::optimize(f = LLSharedBeta,interval=c(0.0001,100),
                            tree=tree, gene.data=gene.data)
  sharedBeta <- sharedBetaFit$minimum

  cat("fit with shared beta =",sharedBeta,"...\n")
  sharedBetaRes <- fitSharedBeta(sharedBeta, tree, gene.data, colSpecies)

  # calculate likelihood ratio test statistic
  LRT <- 2 * (indivBetaRes$ll - sharedBetaRes$ll)


  return( list(indivBetaRes = indivBetaRes,
               sharedBeta = sharedBeta,
               sharedBetaRes = sharedBetaRes,
               LRT = LRT) )
}
