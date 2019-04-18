# This is the main file in which to conduct tests, and calculate the resulting likelihood ratios and chi sqaured values.
# Author: Rori Rohlfs, Lars Gronvold, John Mendoza
# Date: 3/2/2019

library(ape)

source('./scripts/eve-io.R')
source('./scripts/dvdt.R')

divergenceDiversityTest <- function()
{
  # Initialize the tree
  #tree <- read.tree("data/salmonidsPhylo.newick")
  tree <- read.tree("data/examplePhylo2.newick")
  
  # Initialize the gene data
  #gene.data <- getExprData("data/salmonidsBSNsglDupA.tsv")
  gene.data <- getExprData("data/sampleExpr2.dat")
  
  return(betaSharedTest(tree, gene.data))
  
  # Some old debugging code for comparing the results of the betaSharedTest with the dummy data
  #ourIndivBetaResultMatrix <- sapply(indivBetaRes, function(indivBetaResRow) {indivBetaResRow$par} )
  
  #testResMatrix <- t(as.matrix(read.table("data/indivBetaMLParams_trialRun.res", col.names = names(indivBetaRes[[1]]$par))))
  
  #ourIndivBetaResultMatrix / testResMatrix
  
  #ourIndivBetaLLs <- sapply(indivBetaRes, function(indivBetaResRow) {indivBetaResRow$value})
  
  #testResLLs <- read.table("data/indivBetaMLs_trialRun.res")
  #testResLLsVec <- sapply(testResLLs, function(testResLLsRow) {testResLLsRow})
  
  #ourIndivBetaLLs / abs(testResLLs * 1.001)
}


