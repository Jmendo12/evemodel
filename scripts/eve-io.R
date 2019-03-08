# This file contains functions to read input from data files for the EVE-model tests. To include these functions in other files
# code source('eve-io.R') into the file you'd like to use these methods in.
# Author: Rori Rohlfs, Lars Gronvold, John Mendoza
# Data 2/25/19

#Include all our necessary libraries
library(ape)

## Phylogengy Implementation
# Initiliaze the phylogeny data and plot the phylogeny
initializePhylogeny <- function()
{
  # Read in the filename the phylogeny is stored in
  fileName <- paste("./data/", sep = "", readline(prompt = "Enter the file in which the phylogeny data is stored, and include the file extension: "))
  # Create a tree from the file given by the user
  tree <- read.tree(fileName)
  # Reset the margins
  par(mar = c(5, 4, 4, 2) + 0.1)
  # Plot the tree so the edge, node numbers, and tip labels are shown
  plot.phylo(tree, type = "p", show.tip.label = T, label.offset = .05, y.lim = c(.5, 5.5), no.margin = F)
  nodelabels(text = 1 : (Ntip(tree) + Nnode(tree)), node = 1 : (Ntip(tree) + Nnode(tree)), frame = "circle")
  edgelabels(frame = "none", col = "red", adj = c(.5, -.1))
  axisPhylo(1)
  
  return(tree)
}

# Use this function to get the data on the number of samples per species
getIndividuals <- function()
{
  fileName <- paste("./data/", sep = "", readline(prompt = "Enter the file in which the number of samples per species data is stored, and include the file extension: "))
  num.indivs <- scan(fileName, sep = " ")
  return(num.indivs)
}

# Use this function to get the per gene expression data
getExprData <- function(filename)
{
  gene.data <- as.matrix(read.table(filename,header = F,skip = 1,row.names = 1))
  colnames(gene.data) <- scan(filename,nlines = 1,what = character())
 
  return(gene.data)
}