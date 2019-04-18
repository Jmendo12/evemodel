
#### Select a subset of 20 genes from the salmonid data ####

# load salmon example data
source('./scripts/eve-io.R')
tree <- read.tree("data/salmonData/salmonidsPhylo.newick")
gene.data <- getExprData("data/salmonData/salmonidsBSNsglDupA.tsv")

# select the first 20 genes with all values over 0
# (the old implementation may have a problem with negative expression values)
salmon20.data <- gene.data[apply(gene.data>0,1,all),][1:20, ]

# save expression data
write.table(salmon20.data,file = "data/comparison/newFormat/salmonData/salmon20.tsv", quote = F)


#### save in old format ####

# the names of the tips must match the order of the species in the expression data
species <- unique(colnames(salmon20.data))
treeOld <- tree
treeOld$tip.label=match(tree$tip.label,species)
writeLines(c(Ntip(tree),write.tree(treeOld)),"data/comparison/oldEVE/salmonData/salmon.newick")


# Write expression matrix
write(nrow(salmon20.data),file = "data/comparison/oldEVE/salmonData/salmon20.dat")
write.table(salmon20.data,quote = F,sep = " ",col.names = F,file = "data/comparison/oldEVE/salmonData/salmon20.dat",append = T)

# write nindiv
nindiv <- as.vector(table(colnames(salmon20.data))[species])
cat(nindiv, file = "data/comparison/oldEVE/salmonData/salmon.nindiv")
