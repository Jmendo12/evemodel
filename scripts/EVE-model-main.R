# This is the main file in which to conduct tests, and calculate the resulting likelihood ratios and chi sqaured values.
# Author: Rori Rohlfs, Lars Gronvold, John Mendoza
# Date: 3/2/2019

source('./scripts/eve-io.R')
source('./scripts/dvdt.R')
source('./scripts/dvdtsb.R')

divergenceDiversityTest <- function()
{
  # Initialize the tree
  tree <- initializePhylogeny()
  
  # Initialize the number of samples per species
  num.indivs <- getIndividuals()
  
  # Initialize the gene data
  gene.data <- getExprData(num.indivs)
  
  # Calculate the total likelihoods for the null and alternative hypothesis
  ll.IndivBeta <- calculateLLIndivBeta(tree, num.indivs, gene.data)
  ll.SharedBeta <- calculateLLsharedBeta(tree, num.indivs, gene.data)
  
  # Calculate the likelihood ratios - the if statements are to prevent NaN values from being returned from the natural
  # log calculations; these NaNs are returned if the argument passed in to log is negative
  if( ll.IndivBeta < 0 )
  {
    lr <- -2 * log(-ll.IndivBeta) - (-2 * log(ll.SharedBeta))
  }
  else if(ll.SharedBeta < 0)
  {
    lr < - 2 * log(ll.IndivBeta) - (-2 * log(-ll.SharedBeta))
  }
  else
  {
    lr <- -2 * log(ll.IndivBeta) - (-2 * log(ll.SharedBeta))
  }
  
  # If the likelihood ratio is negative, make it positive so that caculating log(lr) does not return NaN
  if(lr < 0)
  {
    lr <- -lr
  }
  
  # Calculate chi sqaured with one degree of freedom
  chi.squared <- 2 * log(lr)
  
  # Since R does not support returning multiple values from a function, save both the likelihood ration value and chi sqaured
  # value to a results file. To view these results in your enviromnent call load("./results/dvdtresults.RData")
  save(lr , chi.squared, file = "./results/dvdtresults.RData")
  
  # Return the chi squared value
  return(chi.squared)
}