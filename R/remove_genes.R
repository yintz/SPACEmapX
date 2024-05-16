
#' @title remove_genes()
#'
#' @description infercnv obj accessor method to remove genes from the matrices
#'
#' @param infercnv_obj infercnv object
#' 
#' @param gene_indices_to_remove matrix indices for genes to remove
#'
#' @return infercnv_obj
#'
#' @keywords internal
#' @noRd
#'

remove_genes <- function(infercnv_obj, gene_indices_to_remove) {

    infercnv_obj@expr.data <- infercnv_obj@expr.data[ -1 * gene_indices_to_remove, , drop=FALSE]
    
    infercnv_obj@count.data <- infercnv_obj@count.data[ -1 * gene_indices_to_remove, , drop=FALSE]

    infercnv_obj@gene_order <- infercnv_obj@gene_order[ -1 * gene_indices_to_remove, , drop=FALSE] 
    infercnv_obj@gene_order[["chr"]] = droplevels(infercnv_obj@gene_order[["chr"]])
    
    validate_infercnv_obj(infercnv_obj)
    
    return(infercnv_obj)
}
