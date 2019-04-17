# Get some expression data from gareths files on the orion server
# The corrsponding phylogeny is in "salmonidsPhylo.newick"
# This data is just for testing.
# TODO: delete this script and data files when proper data has been acquired

library(tidyverse)

load("~/../garethg/SEE/EVE.normalised.exprs.tables.08.11.RData")

# use the between species normalised using single-copy genes as reference (BSNsgl) of dupA data
tbl <- 
  norm.exprs.tables$BSNsgl$dupA %>% 
  drop_na() %>% # Drop rows with NAs
  tibble::column_to_rownames(var = "clan") %>% 
  as.matrix()

colnames(tbl) <- sub("\\..*$","",colnames(tbl))

write.table(tbl,file = "data/salmonidsBSNsglDupA.tsv", quote = F)
