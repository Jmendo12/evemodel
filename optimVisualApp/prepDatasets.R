library(ape)
library(ggtree)
library(tidyverse)

source("./scripts/EVEcore.R")
source('./scripts/eve-io.R')
source("./scripts/twoThetaTest.R")


dataSets <- list(
  salmon20 = list(
    # salmon 20 genes example data
    tree = read.tree("data/comparison/newFormat/salmonData/salmonidsPhylo.newick"),
    gene.data = getExprData("data/comparison/newFormat/salmonData/salmon20.tsv")
  ),
  sim = list(
    # simulated example data
    tree = read.tree("data/comparison/newFormat/simData/examplePhyloNamed.newick"),
    gene.data = getExprData("data/comparison/newFormat/simData/sampleExpr.tsv")    
  )
)

dataSets$salmon20$shiftSpecies <- c("Salp","Okis","Omyk","Ssal")
dataSets$sim$shiftSpecies <- c("D","E")

lapply(dataSets, function(dataSet){
  dataSet$g_phylo <- ggtree(dataSet$tree)
  dataSet$orderedSpcs <- 
    dataSet$g_phylo$data %>% filter(isTip) %>% arrange(y) %>% with(label)
  dataSet$initParams <- initialParamsTwoTheta(dataSet$gene.data,colSpecies = colnames(dataSet$gene.data),shiftSpecies = dataSet$shiftSpecies)
  return(dataSet)
}) -> dataSets

save(dataSets,file = "optimVisualApp/data.RData")
