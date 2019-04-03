library(ape)
source('./scripts/eve-io.R')

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

save(dataSets,file = "optimVisualApp/data.RData")
load("optimVisualApp/data.RData")
