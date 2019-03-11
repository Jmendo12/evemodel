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
  
  # Initialize the number of samples per species
  num.indivs <- as.vector(table(colnames(gene.data)))
  
  # Calculate the total likelihoods for the null and alternative hypothesis
  ll.IndivBeta <- calculateLLIndivBeta(tree, gene.data)
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