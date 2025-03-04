
#' @title apply_median_filtering
#' 
#' @description Apply a median filtering to the expression matrix within each tumor bounds
#'
#' @param SPACEmapX_obj SPACEmapX_object
#' 
#' @param window_size Size of the window side centered on the data point to filter (default = 7).
#' 
#' @param on_observations  boolean (default=TRUE), run on observations data (tumor cells).
#'
#' @param on_references  boolean (default=TRUE), run on references (normal cells).
#'
#' @return SPACEmapX_obj with median filtering applied to observations
#'
#' @export
#'
#' @examples
#' # data(SPACEmapX_data_example)
#' # data(SPACEmapX_annots_example)
#' # data(SPACEmapX_genes_example)
#'
#' # SPACEmapX_object_example <- SPACEmapX::CreateSPACEmapXObject(raw_counts_matrix=SPACEmapX_data_example, 
#' #                                                           gene_order_file=SPACEmapX_genes_example,
#' #                                                           annotations_file=SPACEmapX_annots_example,
#' #                                                           ref_group_names=c("normal"))
#'
#' # SPACEmapX_object_example <- SPACEmapX::run(SPACEmapX_object_example,
#' #                                          cutoff=1,
#' #                                          out_dir=tempfile(), 
#' #                                          cluster_by_groups=TRUE, 
#' #                                          denoise=TRUE,
#' #                                          HMM=FALSE,
#' #                                          num_threads=2,
#' #                                          no_plot=TRUE)
#'
#' data(SPACEmapX_object_example)
#'
#' SPACEmapX_object_example <- SPACEmapX::apply_median_filtering(SPACEmapX_object_example)
#' # plot result object
#'

apply_median_filtering <- function(SPACEmapX_obj,
                                   window_size=7,
                                   on_observations=TRUE,
                                   on_references=TRUE) {

    if (window_size%%2 != 1 | window_size < 2) {
      flog.error("::apply_median_filtering: Error, window_size is an even or < 2. Please specify an odd number >= 3.")
    }
    
    half_window = (window_size - 1) / 2
    
    gene_chr_listing = SPACEmapX_obj@gene_order[[C_CHR]]
    chrs = unlist(unique(gene_chr_listing))
    
    if (on_observations) {
        for (tumor_type in names(SPACEmapX_obj@observation_grouped_cell_indices)) {
            
            tumor_indices_list = SPACEmapX_obj@tumor_subclusters[["subclusters"]][[ tumor_type ]]
            
            for (tumor_indices in tumor_indices_list) {
                for (chr in chrs) {
                    chr_genes_indices = which(gene_chr_listing == chr)
                    working_data = SPACEmapX_obj@expr.data[chr_genes_indices, tumor_indices, drop=FALSE]
                    
                    SPACEmapX_obj@expr.data[chr_genes_indices, tumor_indices] = .median_filter(data=working_data,
                                                                                              window_size=window_size,
                                                                                              half_window=half_window)
                }
            }
        }
    }
    
    if (on_references) {
        for (ref_indices in SPACEmapX_obj@reference_grouped_cell_indices) {
            for (chr in chrs) {
                chr_genes_indices = which(gene_chr_listing == chr)
                working_data = SPACEmapX_obj@expr.data[chr_genes_indices, ref_indices, drop=FALSE]
                
                SPACEmapX_obj@expr.data[chr_genes_indices, ref_indices] = .median_filter(data=working_data,
                                                                                        window_size=window_size,
                                                                                        half_window=half_window)
            }
        }
    }
    
    return(SPACEmapX_obj)
}


.median_filter <- function(data,
                           window_size,
                           half_window) {

    xdim = dim(data)[1]
    ydim = dim(data)[2]
    results = data
    
    # if (xdim >= window_size & ydim >= window_size) {
        for (posx in seq_len(xdim)) {
            posxa <- ifelse(posx <= (half_window + 1), 1, (posx - (half_window + 1)))
            posxb <- ifelse(posx >= (xdim - (half_window + 1)), xdim, (posx + (half_window + 1)))
            for (posy in seq_len(ydim)) {
                posya <- ifelse(posy <= (half_window + 1), 1, (posy - (half_window + 1)))
                posyb <- ifelse(posy >= (ydim - (half_window + 1)), ydim, (posy + (half_window + 1)))
                results[posx, posy] = median(data[posxa:posxb, posya:posyb])
            }
        }
    #}

    return(results)
}
