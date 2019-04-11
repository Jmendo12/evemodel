source('./scripts/dvdt.R')
source('./scripts/eve-io.R')

compareTests <- function()
{
  tree <- read.tree("data/comparison/newFormat/simData/examplePhyloNamed.newick")
  gene.data <- getExprData("data/comparison/newFormat/simData/sampleExpr.tsv")
  
  cat("Timing new implementation \n\n")
  timeNew <- system.time( new_res <- betaSharedTest2(tree, gene.data))
  cat("New implementation took ", timeNew[1:3], " seconds \n\n\n")
  cat("Timing old implementation \n\n")
  timeOld <- system.time( old_res <- betaSharedTest(tree, gene.data))
  cat("Old implementation took ", timeOld[1:3], " seconds \n\n")
  cat("New implementation was ", timeOld[1]/timeNew[1], " times faster than old\n\n")
  
  cat("Return true if number of iterations for each gene with log  transform is less than number of iterations without for individual Beta: \n")
  betterIndiv <- mapply(function(newRow, oldRow) { newRow$counts[1] < oldRow$counts[1] }, 
                        new_res[["indivBetaRes"]], old_res[["indivBetaRes"]])
  cat(betterIndiv, "\n\n")
  
  cat("Return true if number of iterations for each gene with log  transform is less than number of iterations without for shared Beta: \n")
  betterShared <- mapply(function(newRow, oldRow) { newRow$counts[1] < oldRow$counts[1] }, 
                         new_res[["sharedBetaRes"]], old_res[["sharedBetaRes"]])
  cat(betterShared)
}