library(tidyverse)
# get the example data

exprTbl <- readRDS("~/Dropbox/REWIRED project/gitlab/sandve-lab/reEVE/data/BSNormalize/combExprMat.RDS")

colnames(exprTbl) <- c("Drer_1", "Drer_2", "Drer_3", "Drer_4", "Olat_1", "Olat_2", "Olat_3", "Olat_4", "Eluc_1",
                       "Eluc_2", "Eluc_3", "Eluc_4", "Ssal_1", "Ssal_2", "Ssal_3", "Ssal_4", "Salp_1", "Salp_2",
                       "Salp_3", "Salp_4", "Omyk_1", "Omyk_2", "Omyk_3", "Okis_1", "Okis_2", "Okis_3", "Okis_4")
library(tidyverse)

cbind(ID=rownames(exprTbl),exprTbl) %>% 
  as_tibble() %>% 
  na.omit() %>% 
  sample_n(100) %>% 
  write_tsv("exprTbl.tsv")


tree <- read.tree("~/Dropbox/REWIRED project/gitlab/sandve-lab/reEVE/data/from_ortho_pipeline/SpeciesTree_rooted_node_labels.txt")

tree <- keep.tip(tree,c("Drer", "Olat", "Eluc", "Salp", "Okis", "Omyk", "Ssal"))
plot(tree)

tree$node.label <- NULL

write.tree(tree,file = "speciesTree.nwk")
