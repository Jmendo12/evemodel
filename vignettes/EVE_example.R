## ----setup, include=FALSE------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

library(evemodel)

## ------------------------------------------------------------------------
# read tab separated table
exprTbl <- read.delim(system.file("extdata", "exprTbl.tsv", package = "evemodel"))

head(exprTbl)

## ------------------------------------------------------------------------
# The first column in the table is the ortholog group ID
# which will be the rownames of the matrix.
# TODO: The rownames are currently not returned in the results, but they should.
exprMat <- as.matrix(exprTbl[,-1])
rownames(exprMat) <- exprTbl$ID

## ------------------------------------------------------------------------
library(ape)

speciesTree <- read.tree(system.file("extdata", "speciesTree.nwk", package = "evemodel"))

# plot the species tree
plot(speciesTree)
add.scale.bar(x=0,y=7,length = 0.1)

## ------------------------------------------------------------------------
# the species names in the tree is given by the tip.labels
speciesTree$tip.label

# the columns in the expression matrix are:
colnames(exprMat)

# remove the trailing number so that we get a vector with the species for each column
colSpecies <- sub("_.*$","",colnames(exprMat))

colSpecies

## ------------------------------------------------------------------------
# plot likelihood ratio test statistic histogram
hist(res$LRT,freq = F)

# Plot the chi-squared distribution with one degree of freedom
x = seq(0.5,10,length.out = 100)
y = dchisq(x,df = 1)
lines(x,y,col="red")

## ------------------------------------------------------------------------
pval = pchisq(res$LRT,df = 1,lower.tail = F)

## ------------------------------------------------------------------------
res$sharedBeta

## ------------------------------------------------------------------------
head(res$sharedBetaRes$par)

## ------------------------------------------------------------------------
head(res$indivBetaRes$par)

## ------------------------------------------------------------------------
# get median of each parameter
apply(res$indivBetaRes$par,2,median)

## ------------------------------------------------------------------------
set.seed(123) # we want the same simulation each time...

simData <- 
  simOneTheta(n = 100, tree = speciesTree, colSpecies = colSpecies,
              theta = 2, sigma2 = 100, alpha = 25, beta = 0.3)

## ------------------------------------------------------------------------
# plot likelihood ratio test statistic histogram
hist(resSim$LRT,freq = F)

# Plot the chi-squared distribution with one degree of freedom
x = seq(0.2,10,length.out = 100)
y = dchisq(x,df = 1)
lines(x,y,col="red")

