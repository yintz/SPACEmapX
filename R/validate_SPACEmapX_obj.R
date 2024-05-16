
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

validate_SPACEmapX_obj <- function(SPACEmapX_obj) {

    flog.info("validating SPACEmapX_obj")

    if (isTRUE(all.equal(rownames(SPACEmapX_obj@expr.data), rownames(SPACEmapX_obj@gene_order)))) {
            # all good.
                    return();

    }
        else {

        flog.error("hmm.... rownames(SPACEmapX@expr.data != rownames(SPACEmapX@gene_order))")
                broken.SPACEmapX_obj = SPACEmapX_obj
                        save('broken.SPACEmapX_obj', file="broken.SPACEmapX_obj")

    }


    genes = setdiff(rownames(SPACEmapX_obj@expr.data), rownames(SPACEmapX_obj@gene_order))
        if (length(genes) != 0) {
                flog.error(paste("The following genes are in SPACEmapX_obj@expr.data and not @gene_order:", paste(genes, collapse=","),
                                         sep=" "))

    }

    genes = setdiff(rownames(SPACEmapX_obj@gene_order), rownames(SPACEmapX_obj@expr.data))
        if (length(genes) != 0) {
                flog.error(paste("The following genes are in @gene_order and not SPACEmapX_obj@expr.data:", paste(genes, collapse=","),
                                         sep=" "))

    }

    stop("Problem detected w/ SPACEmapX_obj")

}


get_cell_name_by_grouping <- function(SPACEmapX_obj) {

    cell_name_groupings = list()
    
    groupings = c(SPACEmapX_obj@reference_grouped_cell_indices, SPACEmapX_obj@observation_grouped_cell_indices)

    for (group_name in names(groupings)) {

        cell_names = colnames(SPACEmapX_obj@expr.data[, groupings[[ group_name ]] ] )

        cell_name_groupings[[ group_name ]] = cell_names

    }

    return(cell_name_groupings)
}


has_reference_cells <- function(SPACEmapX_obj) {
    return(length(SPACEmapX_obj@reference_grouped_cell_indices) != 0)
}
