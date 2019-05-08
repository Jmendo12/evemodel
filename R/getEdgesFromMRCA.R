#' Get the edges for a clade in a phylogeny given the tips
#' 
#'
#' @param tree Phylogeny
#' @param tips Name of the tips that defines the clade
#' @param includeEdgeToMRCA If TRUE, also include the edge leading to the root of the clade
#'
#' @return vector of edge indices
#' @export
getEdgesFromMRCA <- function(tree, tips, includeEdgeToMRCA){
  mrcaNode <- ape::getMRCA(tree, tips)
  
  getDescendantEdgesRecursive <- function(nodes){
    descendantEdges <- which(tree$edge[ ,1] %in% nodes)
    if( length(descendantEdges) > 0){
      descendantEdges <- c(descendantEdges, Recall(nodes = tree$edge[descendantEdges,2]))
    }
    return(descendantEdges)
  }
  
  edges <- getDescendantEdgesRecursive(mrcaNode)
  
  if( includeEdgeToMRCA ){
    edges <- c(edges,match(mrcaNode, tree$edge[,2]))
  } 
  
  return( edges )
}