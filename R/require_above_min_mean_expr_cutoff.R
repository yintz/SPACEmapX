
#' @title require_above_min_mean_expr_cutoff ()
#'
#' @description Filters out genes that have fewer than the corresponding mean value across all cell values.
#'
#' @param infercnv_obj  infercnv_object
#'
#' @param min_mean_expr_cutoff  the minimum mean value allowed for a gene to be retained in the expression matrix.
#'
#' @return infercnv_obj  the infercnv_object with lowly or unexpressed genes removed.
#'
#' @keywords internal
#' @noRd
#'






require_above_min_mean_expr_cutoff <- function(infercnv_obj, min_mean_expr_cutoff) {
    
    flog.info(paste("::above_min_mean_expr_cutoff:Start", sep=""))
    
    
    indices <-.below_min_mean_expr_cutoff(infercnv_obj@expr.data, min_mean_expr_cutoff)
    if (length(indices) > 0) {
        flog.info(sprintf("Removing %d genes from matrix as below mean expr threshold: %g",
                          length(indices), min_mean_expr_cutoff))
        
        infercnv_obj <- remove_genes(infercnv_obj, indices)
        
        expr_dim = dim(infercnv_obj@expr.data)
        flog.info(sprintf("There are %d genes and %d cells remaining in the expr matrix.",
                          expr_dim[1], expr_dim[2]))
        
    }
    
    return(infercnv_obj)
    
}
