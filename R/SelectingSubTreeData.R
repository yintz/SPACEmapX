#' Selecting Subtree Data for Twig Selection: this selects a number of barcoded spots from a inferCNV dendrogram object for further analysis.
#'
#' SelectingSubTreeData()
#' 
#' @param SubtreeObject A dendrogram, phylo object created by subtrees(as.phylo([dendogram.txt]))
#' @param TwigOfInterest A numerical integer corresponding to a phylogram/dendogram twig of interest
#' @return A specific subtree twig
#' @examples
#' SelectingSubTreeData(my.subtrees, 4617)

SelectingSubTreeData <- function(SubtreeObject, TwigOfInterest) {
  tree_twig <- SubtreeObject[[TwigOfInterest]]
  output <- tree_twig$tip.label
  output <- as.data.frame(output)
  output <- output %>%
    mutate(Twig = paste0("Twig_", TwigOfInterest))
  names(output)[1] <- "Barcode"
  return(output)
}
