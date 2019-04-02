# This file contains functions to read input from data files for the EVE-model tests. To include these functions in other files
# code source('eve-io.R') into the file you'd like to use these methods in.
# Author: Rori Rohlfs, Lars Gronvold, John Mendoza
# Data 2/25/19


# Use this function to get the per gene expression data
getExprData <- function(filename)
{
  gene.data <- as.matrix(read.table(filename,header = F,skip = 1,row.names = 1))
  colnames(gene.data) <- scan(filename,nlines = 1,what = character())
  
  # Filter out the data that causes errors and warnings to arise from the optimization
  return(gene.data[apply(gene.data, 1, var) > 0, ])
}

getSpeciesToShift <- function()
{
  shiftedSpecies <- readline(prompt = "Enter the species names to shift: ")
}