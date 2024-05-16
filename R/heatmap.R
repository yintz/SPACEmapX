#!/usr/bin/env Rscript


# Returns the color palette for contigs.
#
# Returns:
# Color Palette
get_group_color_palette <- function(){
    return(colorRampPalette(RColorBrewer::brewer.pal(12,"Set3")))
}


#' @description Formats the data and sends it for plotting.
#'
#' @title Plot the matrix as a heatmap, with cells as rows and genes as columns, ordered according to chromosome
#'
#' @param infercnv_obj infercnv object
#' @param out_dir Directory in which to save pdf and other output.
#' @param title Plot title.
#' @param obs_title Title for the observations matrix.
#' @param ref_title Title for the reference matrix.
#' @param cluster_by_groups Whether to cluster observations by their annotations or not. Using this ignores k_obs_groups.
#' @param cluster_references Whether to cluster references within their annotations or not. (dendrogram not displayed)
#' @param plot_chr_scale Whether to scale the chromosme width on the heatmap based on their actual size rather than just the number of expressed genes.
#' @param chr_lengths A named list of chromsomes lengths to use when plot_chr_scale=TRUE, or else chromosome size is assumed to be the last chromosome's stop position + 10k bp
#' @param k_obs_groups Number of groups to break observation into.
#' @param contig_cex Contig text size. 
#' @param x.center Value on which to center expression.
#' @param x.range vector containing the extreme values in the heatmap (ie. c(-3,4) )
#' @param hclust_method Clustering method to use for hclust.
#' @param custom_color_pal Specify a custom set of colors for the heatmap. 
#'                         Has to be in the shape color.palette(c("darkblue", "white", "darkred"),
#'                                                              c(2, 2))
#' @param color_safe_pal Logical indication of using a color blindness safe palette.
#' @param output_filename Filename to save the figure to.
#' @param output_format format for heatmap image file (default: 'png'), options('png', 'pdf', NA)
#'                      If set to NA, will print graphics natively
#' @param png_res Resolution for png output.
#' @param dynamic_resize Factor (>= 0) by which to scale the dynamic resize of the observation 
#'                       heatmap and the overall plot based on how many cells there are.
#'                       Default is 0, which disables the scaling. Try 1 first if you want to enable.
#' @param ref_contig If given, will focus cluster on only genes in this contig.
#' @param write_expr_matrix Includes writing a matrix file containing the expression data that is plotted in the heatmap.
#' @param write_phylo Write newick strings of the dendrograms displayed on the left side of the heatmap to file.
#' @param useRaster Whether to use rasterization for drawing heatmap. Only disable if it produces an error as it is much faster than not using it.
#' 
#' @return A list of all relevent settings used for the plotting to be able to reuse them in another plot call while keeping consistant plotting settings, most importantly x.range.
#'
#'
#' @export
#'
#' @examples
#' # data(infercnv_data_example)
#' # data(infercnv_annots_example)
#' # data(infercnv_genes_example)
#'
#' # infercnv_object_example <- infercnv::CreateInfercnvObject(raw_counts_matrix=infercnv_data_example, 
#' #                                                           gene_order_file=infercnv_genes_example,
#' #                                                           annotations_file=infercnv_annots_example,
#' #                                                           ref_group_names=c("normal"))
#'
#' # infercnv_object_example <- infercnv::run(infercnv_object_example,
#' #                                          cutoff=1,
#' #                                          out_dir=tempfile(), 
#' #                                          cluster_by_groups=TRUE, 
#' #                                          denoise=TRUE,
#' #                                          HMM=FALSE,
#' #                                          num_threads=2,
#' #                                          no_plot=TRUE)
#'
#' data(infercnv_object_example)
#'
#' plot_cnv(infercnv_object_example,
#'          out_dir=tempfile(),
#'          obs_title="Observations (Cells)",
#'          ref_title="References (Cells)",
#'          cluster_by_groups=TRUE,
#'          x.center=1,
#'          x.range="auto",
#'          hclust_method='ward.D',
#'          color_safe_pal=FALSE,
#'          output_filename="infercnv",
#'          output_format="png",
#'          png_res=300,
#'          dynamic_resize=0
#'          )
#'



C_CHR <- "chr"
C_START <- "start"
C_STOP <- "stop"
C_HCLUST_METHODS <- c("ward.D", "ward.D2", "single", "complete", "average", "mcquitty", "median", "centroid")
C_OUTPUT_FORMAT <- c("pdf", "png")




plot_cnv <- function(infercnv_obj,
                     out_dir=".",
                     title="inferCNV",
                     obs_title="Observations (Cells)",
                     ref_title="References (Cells)",
                     cluster_by_groups=TRUE,
                     cluster_references=TRUE,
                     plot_chr_scale=FALSE,
                     chr_lengths=NULL,
                     k_obs_groups = 1,
                     contig_cex=1,
                     x.center=mean(infercnv_obj@expr.data),
                     x.range="auto", #NA,
                     hclust_method='ward.D',
                     custom_color_pal=NULL,
                     color_safe_pal=FALSE,
                     output_filename="infercnv",
                     output_format="png", #pdf, png, NA
                     png_res=300,
                     dynamic_resize=0,
                     ref_contig = NULL,
                     write_expr_matrix=FALSE,
                     write_phylo=FALSE,
                     useRaster=TRUE) {


    # arg validations
    if (! hclust_method %in% C_HCLUST_METHODS) {
        stop(sprintf("Error, hclust_method: %s is not supported", hclust_method))
    }
    if ( (! is.na(output_format) ) & (!  output_format %in% C_OUTPUT_FORMAT) )  {
        stop(sprintf("Error, output_format: %s is not supported", output_format) )
    }
    
    if(!file.exists(out_dir)){
        dir.create(out_dir)
    }
    
    plot_data = infercnv_obj@expr.data
    
    flog.info(paste("::plot_cnv:Start", sep=""))
    flog.info(paste("::plot_cnv:Current data dimensions (r,c)=",
                           paste(dim(plot_data), collapse=","),
                           " Total=", sum(plot_data, na.rm=TRUE),
                           " Min=", min(plot_data, na.rm=TRUE),
                           " Max=", max(plot_data, na.rm=TRUE),
                           ".", sep=""))
    flog.info(paste("::plot_cnv:Depending on the size of the matrix",
                           " this may take a moment.",
                           sep=""))


    
    if (write_expr_matrix) {
        expr_dat_file <- paste(out_dir, paste("expr.", output_filename, ".dat", sep=""), sep="/")

        if ("matrix" %in% is(plot_data)) {
            write.table(as.matrix(plot_data), file=expr_dat_file, quote=FALSE, sep="\t")
        }
        
    }
    
    if (! any(is.na(x.range))) {

        if ( (length(x.range) == 1) & (x.range[1] == "auto") ) {

            # examine distribution of data that's off-center, since much of the center could
            # correspond to a mass of data that has been wiped out during noise reduction
            quantiles = quantile(plot_data[plot_data != x.center], c(0.01, 0.99))

            # determine max distance from the center.
            delta = max( abs( c(x.center - quantiles[1],  quantiles[2] - x.center) ) )
            low_threshold = x.center - delta
            high_threshold = x.center + delta
            x.range = c(low_threshold, high_threshold)
            
            flog.info(sprintf("plot_cnv(): auto thresholding at: (%f , %f)", low_threshold, high_threshold))
            
        } else {
        
            # use defined values
            low_threshold = x.range[1]
            high_threshold = x.range[2]
            
            if (low_threshold > x.center | high_threshold < x.center | low_threshold >= high_threshold) {
                stop(paste("Error, problem with relative values of x.range: ", x.range, ", and x.center: ", x.center))
            }
        }

        plot_data[plot_data < low_threshold] <- low_threshold
        plot_data[plot_data > high_threshold] <- high_threshold
        
        infercnv_obj@expr.data <- plot_data  #because used again below...
        
    }
    
    
    # Contigs
    contigs = infercnv_obj@gene_order[[C_CHR]]
    unique_contigs <- unique(contigs)
    n_contig <- length(unique_contigs)
    ct.colors <- get_group_color_palette()(n_contig)
    names(ct.colors) <- unique_contigs

    # Select color palette
    if (!is.null(custom_color_pal)) {
        custom_pal = custom_color_pal
    } else if (color_safe_pal == FALSE){
        custom_pal <- color.palette(c("darkblue", "white", "darkred"),
                                    c(2, 2))
    } else {
        custom_pal <- color.palette(c("purple3", "white", "darkorange2"),
                                    c(2, 2))
    }

    
    ## Row separation based on reference
    ref_idx <- NULL
    if (has_reference_cells(infercnv_obj)) {
        ref_idx <- unlist(infercnv_obj@reference_grouped_cell_indices)
        ref_idx = ref_idx[order(ref_idx)]
    }
    
    # Column seperation by contig and label axes with only one instance of name
    contig_tbl <- table(contigs)[unique_contigs]
    col_sep <- cumsum(contig_tbl)
    col_sep <- col_sep[-1 * length(col_sep)]   ## FIXME:  removing last entry?
    # These labels are axes labels, indicating contigs at the first column only
    # and leaving the rest blank.
    # contig_labels <- c()
    # contig_names <-c()
    # for (contig_name in names(contig_tbl)){
    #     contig_labels <- c(contig_labels,
    #                        contig_name,
    #                        rep("", contig_tbl[contig_name] - 1))
    #     contig_names <- c(contig_names,rep(contig_name,contig_tbl[contig_name]))
    # }
    contig_labels = names(contig_tbl)
    contig_names = unlist(lapply(contig_labels, function(contig_name) {rep(contig_name, contig_tbl[contig_name])}))

    # contig_sizes = vapply(contig_labels, function(contig_name) {contig_tbl[contig_name]}, integer(1))
    # contig_positions =  c(1, cumsum(contig_sizes[seq_len(length(contig_sizes) - 1)]) + 1)

    # Calculate how many rows will be made for the number of columns in the grouping key
    grouping_key_coln <- c()
    obs_annotations_names <- names(infercnv_obj@observation_grouped_cell_indices)


   
    # obs_annotations_groups: integer vec named by cells, set to index according to category name vec above.
    obs_annotations_groups = rep(-1, length(colnames(infercnv_obj@expr.data))) # init
    names(obs_annotations_groups) = colnames(infercnv_obj@expr.data)
    obs_index_groupings = infercnv_obj@observation_grouped_cell_indices
    counter <- 1
    for (obs_index_group in obs_index_groupings) {
        obs_annotations_groups[ obs_index_group ] <- counter
        counter <- counter + 1
    }
    # restrict to just the obs indices

    if (! is.null(ref_idx)) {
        obs_annotations_groups <- obs_annotations_groups[-ref_idx]
    }
    
    if (is.null(dynamic_resize) | dynamic_resize < 0) {
        flog.warn(paste("invalid dynamic_resize value: ", dynamic_resize, sep=""))
        dynamic_resize = 0
    }
    dynamic_extension = 0
    nobs = length(unlist(infercnv_obj@observation_grouped_cell_indices))
    if (nobs > 200) {
        dynamic_extension = dynamic_resize * 3.6 * (nobs - 200)/200 
    }

    grouping_key_coln[1] <- floor(123/(max(nchar(obs_annotations_names)) + 6))  ## 123 is the max width in number of characters, 4 is the space taken by the color box itself and the spacing around it
    if (grouping_key_coln[1] < 1) {
        grouping_key_coln[1] <- 1
    }

    name_ref_groups = names(infercnv_obj@reference_grouped_cell_indices)
    if (is.null(name_ref_groups)) {
        grouping_key_coln[2] = 1
    } else {
        grouping_key_coln[2] <- floor(123/(max(nchar(name_ref_groups)) + 6))  ## 123 is the max width in number of characters, 4 is the space taken by the color box itself and the spacing around it
        if (grouping_key_coln[2] < 1) {
            grouping_key_coln[2] <- 1
        }
    }

    grouping_key_rown <- c()
    grouping_key_rown[1] <- ceiling(length(obs_annotations_names)/grouping_key_coln[1])
    grouping_key_rown[2] <- ceiling(length(name_ref_groups)/grouping_key_coln[2])
    # Calculate how much bigger the output needs to be to accodomate for the grouping key
    grouping_key_height <- c((grouping_key_rown[2] + 2) * 0.175, (grouping_key_rown[1] + 3) * 0.175)

    # Rows observations, Columns CHR
    if (! is.na(output_format)) {
        if (output_format == "pdf") {
            pdf(paste(out_dir, paste(output_filename, ".pdf", sep=""), sep="/"),
                useDingbats=FALSE,
                width=11,
                height=(8.22 + sum(grouping_key_height)) + dynamic_extension,
                paper="special")
        } else if (output_format == "png") {
            png_height = 8.22 + sum(grouping_key_height) + dynamic_extension
            if (!is.null(getOption("bitmapType")) && (getOption("bitmapType") == "cairo") & (png_height > 32768/png_res)) {  # 32768 is the pixel limit for cairo backend
                png_height = round((32767/png_res) - 5*10^(-3), 2) # floor() with 2 decimals
                dynamic_extension = png_height - 8.22 - sum(grouping_key_height)
                flog.warn(paste0("Requested PNG output height too big at the current resolution, ",
                    "using the max height instead. (cairo seems to have a size limit of 32767 (2^15-1) pixels ",
                    "per dimension and 49151 (2^15+2^14-1)pixels for the sum of dimensions)"))
            }
            png(paste(out_dir, paste(output_filename, ".png", sep=""), sep="/"),
                width=10,
                height=png_height,
                units="in",
                res=png_res)
        }
    }
    


    # Plot observations
    ## Make Observation Samples
    ## Remove observation col names, too many to plot
    ## Will try and keep the reference names
    ## They are more informative anyway
    obs_data <- infercnv_obj@expr.data
    if (!is.null(ref_idx)){
        obs_data <- plot_data[, -ref_idx, drop=FALSE]
        if (ncol(obs_data) == 1) {
            # hack for dealing with single entries
            plot_data <- cbind(obs_data, obs_data)
            names(obs_data) <- c("", names(obs_data)[1])
        }
    }
    
    obs_data <- t(obs_data)

    # Subsample the data to only the references and update the ref_group indexes to match their new indexes
    # ref_data_t <- plot_data[, ref_idx, drop=FALSE]
    ref_data_t <- NULL
    updated_ref_groups <- list()
    current_ref_count <- 1
    current_grp_idx <- 1
    plot_data <-infercnv_obj@expr.data
    ref_groups = infercnv_obj@reference_grouped_cell_indices
    for (ref_grp in ref_groups) {
        ref_data_t <- cbind(ref_data_t, plot_data[, ref_grp, drop=FALSE])
        updated_ref_groups[[current_grp_idx]] = seq(current_ref_count, current_ref_count + length(ref_grp) - 1)
        current_ref_count <- current_ref_count + length(ref_grp)
        current_grp_idx <- current_grp_idx + 1
    }
    ref_groups <- updated_ref_groups
    
    nb_breaks <- 16
    # breaksList_t <-
    #     seq(min(min(obs_data, na.rm=TRUE), min(ref_data_t, na.rm=TRUE)),
    #     max(max(obs_data,na.rm=TRUE), max(ref_data_t, na.rm=TRUE)),
    #     length.out=nb_breaks)
    breaksList_t <- seq(x.range[1], x.range[2], length.out=nb_breaks)


    gene_position_breaks = NULL
    if (plot_chr_scale) {
        # gene table to heatmap width
        chr_name_list = unique(infercnv_obj@gene_order[["chr"]])
        
        # get average distance for tail end? 
        # optionally give vector of chr lengths
        # if chr_lengths not given, take furthest gene and add 10k bp
        if (is.null(chr_lengths)) {
            chr_lengths = c()
            for (chr_name in chr_name_list) {
                chr_lengths = c(chr_lengths, max(infercnv_obj@gene_order$stop[which(infercnv_obj@gene_order$chr == chr_name)]) + 10000)
            }
            names(chr_lengths) = chr_name_list
        }
        gene_position_breaks = vector(mode="integer", length=(length(unlist(infercnv_obj@gene_order$chr)) + 1))

        sum_previous_contigs = 0
        gene_position_breaks[1] = 1
        current_idx = 2
        col_sep_idx = 1

        for (chr_name in chr_name_list) {
            index_pos = which(infercnv_obj@gene_order$chr == chr_name)
            latest_position = 1
            if (length(index_pos) > 1) {
                for (i in index_pos[2:length(index_pos)]) {
                    gene_position_breaks[current_idx] = sum_previous_contigs + ((latest_position + infercnv_obj@gene_order$start[i]) / 2)
                    latest_position = max(infercnv_obj@gene_order$stop[i], latest_position)
                    current_idx = current_idx + 1
                }
            }
            gene_position_breaks[current_idx] = sum_previous_contigs + chr_lengths[chr_name]
            current_idx = current_idx + 1
            sum_previous_contigs = sum_previous_contigs + chr_lengths[chr_name]


            if (col_sep_idx != length(chr_name_list)) {
                col_sep[col_sep_idx] = sum_previous_contigs
                col_sep_idx = col_sep_idx + 1
            }
        }

        overlap_pos = which(gene_position_breaks[1:(length(gene_position_breaks)-1)] == gene_position_breaks[2:length(gene_position_breaks)])
        if (any(overlap_pos)) {
            gene_position_breaks[overlap_pos+1] = gene_position_breaks[overlap_pos] + 1
        }
    }



    # Create file base for plotting output
    force_layout <- .plot_observations_layout(grouping_key_height=grouping_key_height, dynamic_extension=dynamic_extension)
    .plot_cnv_observations(infercnv_obj=infercnv_obj,
                          obs_data=obs_data,
                          file_base_name=out_dir,
                          do_plot=!is.na(output_format),
                          write_expr_matrix=write_expr_matrix,
                          write_phylo=write_phylo,
                          output_filename_prefix=output_filename,
                          cluster_contig=ref_contig,
                          contigs=contigs,
                          contig_colors=ct.colors[contigs],
                          contig_labels=contig_labels,
                          contig_names=contig_names,
                          col_pal=custom_pal,
                          contig_seps=col_sep,
                          num_obs_groups=k_obs_groups,
                          obs_annotations_groups=obs_annotations_groups,
                          obs_annotations_names=obs_annotations_names,
                          grouping_key_coln=grouping_key_coln[1],
                          cluster_by_groups=cluster_by_groups,
                          cnv_title=title,
                          cnv_obs_title=obs_title,
                          contig_lab_size=contig_cex,
                          breaksList=breaksList_t,
                          gene_position_breaks=gene_position_breaks,
                          x.center=x.center,
                          hclust_method=hclust_method,
                          layout_lmat=force_layout[["lmat"]],
                          layout_lhei=force_layout[["lhei"]],
                          layout_lwid=force_layout[["lwid"]],
                          useRaster=useRaster)
    obs_data <- NULL

    if(!is.null(ref_idx)){
        .plot_cnv_references(infercnv_obj=infercnv_obj,
                            ref_data=ref_data_t,
                            ref_groups=ref_groups,
                            name_ref_groups=name_ref_groups,
                            cluster_references=cluster_references,
                            hclust_method=hclust_method,
                            grouping_key_coln=grouping_key_coln[2],
                            col_pal=custom_pal,
                            contig_seps=col_sep,
                            file_base_name=out_dir,
                            do_plot=!is.na(output_format),
                            write_expr_matrix=write_expr_matrix,
                            output_filename_prefix=output_filename,
                            cnv_ref_title=ref_title,
                            breaksList=breaksList_t,
                            gene_position_breaks=gene_position_breaks,
                            x.center=x.center,
                            layout_add=TRUE,
                            useRaster=useRaster)
    }
    if (! is.na(output_format)) {
        dev.off()
    }

    return(list(cluster_by_groups = cluster_by_groups,
               k_obs_groups = k_obs_groups,
               contig_cex = contig_cex,
               x.center = x.center,
               x.range = x.range,
               hclust_method = hclust_method,
               color_safe_pal = color_safe_pal,
               output_format = output_format,
               png_res = png_res,
               dynamic_resize = dynamic_resize))
}

















