
#' @title validate_infercnv_obj()
#'
#' @description validate an infercnv_obj
#' ensures that order of genes in the @gene_order slot match up perfectly with the gene rows in the @expr.data matrix.
#' Otherwise, throws an error and stops execution.
#'
#' @param infercnv_obj infercnv_object
#'
#' @return none
#'

validate_infercnv_obj <- function(infercnv_obj) {

    flog.info("validating infercnv_obj")

    if (isTRUE(all.equal(rownames(infercnv_obj@expr.data), rownames(infercnv_obj@gene_order)))) {
            # all good.
                    return();

    }
        else {

        flog.error("hmm.... rownames(infercnv_obj@expr.data != rownames(infercnv_obj@gene_order))")
                broken.infercnv_obj = infercnv_obj
                        save('broken.infercnv_obj', file="broken.infercnv_obj")

    }


    genes = setdiff(rownames(infercnv_obj@expr.data), rownames(infercnv_obj@gene_order))
        if (length(genes) != 0) {
                flog.error(paste("The following genes are in infercnv_obj@expr.data and not @gene_order:", paste(genes, collapse=","),
                                         sep=" "))

    }

    genes = setdiff(rownames(infercnv_obj@gene_order), rownames(infercnv_obj@expr.data))
        if (length(genes) != 0) {
                flog.error(paste("The following genes are in @gene_order and not infercnv_obj@expr.data:", paste(genes, collapse=","),
                                         sep=" "))

    }

    stop("Problem detected w/ infercnv_obj")

}


get_cell_name_by_grouping <- function(infercnv_obj) {

    cell_name_groupings = list()
    
    groupings = c(infercnv_obj@reference_grouped_cell_indices, infercnv_obj@observation_grouped_cell_indices)

    for (group_name in names(groupings)) {

        cell_names = colnames(infercnv_obj@expr.data[, groupings[[ group_name ]] ] )

        cell_name_groupings[[ group_name ]] = cell_names

    }

    return(cell_name_groupings)
}


has_reference_cells <- function(infercnv_obj) {
    return(length(infercnv_obj@reference_grouped_cell_indices) != 0)
}
