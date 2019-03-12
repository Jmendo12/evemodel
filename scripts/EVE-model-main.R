# This is the main file in which to conduct tests, and calculate the resulting likelihood ratios and chi sqaured values.
# Author: Rori Rohlfs, Lars Gronvold, John Mendoza
# Date: 3/2/2019

library(ape)

source('./scripts/eve-io.R')
source('./scripts/dvdt.R')
source('./scripts/dvdtsb.R')

divergenceDiversityTest <- function()
{
  # Initialize the tree
  tree <- read.tree("data/examplePhylo2.newick")
  
  # Initialize the gene data
  gene.data <- getExprData("data/sampleExpr2.dat")
  
  # Calculate the total likelihoods for the null and alternative hypothesis
  indivBetaRes <- calculateLLIndivBeta(tree, gene.data)
  
  ourIndivBetaResultMatrix <- sapply(indivBetaRes, function(indivBetaResRow) {indivBetaResRow$par} )
  
  testResMatrix <- t(as.matrix(read.table("data/indivBetaMLParams_trialRun.res", col.names = names(indivBetaRes[[1]]$par))))
  
  ourParamResultMatrix <= testResMatrix * 1.001
  
  ourIndivBetaLLs <- sapply(indivBetaRes, function(indivBetaResRow) {indivBetaResRow$value})
  
  testResLLs <- read.table("data/indivBetaMLs_trialRun.res")
  testResLLsVec <- sapply(testResLLs, function(testResLLsRow) {testResLLsRow})
  
  ourIndivBetaLLs <= abs(testResLLs * 1.001)
  
  ll.SharedBeta <- calculateLLsharedBeta(tree, num.indivs, gene.data)
  
  # Calculate the likelihood ratios
  likelihood.ratio <- abs(-2 * log(ll.IndivBeta) - (-2 * log(ll.SharedBeta)))
  
  # Calculate chi sqaured with one degree of freedom
  chi.squared <- 2 * log(likelihood.ratio)
  
  # Since R does not support returning multiple values from a function, save both the likelihood ration value and chi sqaured
  # value to a results file. To view these results in your enviromnent call load("./results/dvdtresults.RData")
  save(likelihood.ratio , chi.squared, file = "./results/dvdtresults.RData")
  
  # Return the chi squared value
  return(chi.squared)
}