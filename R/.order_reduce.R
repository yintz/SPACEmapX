
# Order the data and subset the data to data in the genomic position file.
#
# Args:
# @param data Data (expression) matrix where the row names should be in
#                 the row names of the genomic_position file.
# @param genomic_position Data frame read in from the genomic position file
#
# @return Returns a matrix of expression in the order of the
#            genomic_position file. NULL is returned if the genes in both
#            data parameters do not match.
#





.order_reduce <- function(data, genomic_position){
    flog.info(paste("::order_reduce:Start.", sep=""))
    ret_results <- list(expr=NULL, order=NULL, chr_order=NULL)
    if (is.null(data) || is.null(genomic_position)){
        return(ret_results)
    }

    # Drop pos_gen entries that are position 0
    remove_by_position <- -1 * which(genomic_position[2] + genomic_position[3] == 0)
    if (length(remove_by_position)) {
        flog.debug(paste("::process_data:order_reduce: removing genes specified by pos == 0, count: ",
                                length(remove_by_position), sep=""))

        genomic_position <- genomic_position[remove_by_position, , drop=FALSE]
    }

    # Reduce to genes in pos file

    flog.debug(paste("::process_data:order_reduce: gene identifers in expression matrix: ",
                            row.names(data), collapse="\n", sep=""))
    flog.debug(paste("::process_data:order_reduce: gene identifers in genomic position table: ",
                            row.names(data), collapse="\n", sep=""))



    keep_genes <- intersect(row.names(data), row.names(genomic_position))
        
    flog.debug(paste("::process_data:order_reduce: keep_genes size: ", length(keep_genes),
                     sep=""))
    
    # Keep genes found in position file
    if (length(keep_genes)) {
        ret_results$expr <- data[match(keep_genes, rownames(data)), , drop=FALSE]
        ret_results$order <- genomic_position[match(keep_genes, rownames(genomic_position)), , drop=FALSE]
    } else {
        flog.info(paste("::process_data:order_reduce:The position file ",
                        "and the expression file row (gene) names do not match."))
        return(list(expr=NULL, order=NULL, chr_order=NULL))
    }
    
    ## ensure expr and order match up!
    if (isTRUE(all.equal(rownames(ret_results$expr), rownames(ret_results$order)))) {
        flog.info(".order_reduce(): expr and order match.")
    }
    else {
        stop("Error, .order_reduce(): expr and order don't match! must debug")
    }
        
    # Set the chr to factor so the order can be arbitrarily set and sorted.
    chr_levels <- unique(genomic_position[[C_CHR]])
    ret_results$order[[C_CHR]] <- factor(ret_results$order[[C_CHR]],
                                   levels=chr_levels)
    
    # Sort genomic position file and expression file to genomic position file
    # Order genes by genomic region
    order_names <- row.names(ret_results$order)[with(ret_results$order, order(chr,start,stop))]
    
    ret_results$expr <- ret_results$expr[match(order_names, rownames(ret_results$expr)), , drop=FALSE] #na.omit is to rid teh duplicate gene entries (ie. Y_RNA, snoU13, ...) if they should exist.
    
    # This is the contig order, will be used in visualization.
    # Get the contig order in the same order as the genes.
    ret_results$order <- ret_results$order[match(order_names, rownames(ret_results$order)), , drop=FALSE]
    ret_results$chr_order <- ret_results$order[1]
    
    # Remove any gene without position information
    # Genes may be sorted correctly by not have position information
    # Here they are removed.
    flog.info(paste("::process_data:order_reduce:Reduction from positional ",
                           "data, new dimensions (r,c) = ",
                           paste(dim(data), collapse=","),
                           " Total=", sum(data),
                           " Min=", min(data),
                           " Max=", max(data),
                           ".", sep=""))
    flog.debug(paste("::process_data:order_reduce end."))
    return(ret_results)
}

