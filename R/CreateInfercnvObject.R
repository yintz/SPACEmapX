
#' @title CreateInfercnvObject
#'
#' @param raw_counts_matrix  the matrix of genes (rows) vs. cells (columns) containing the raw counts
#'                           If a filename is given, it'll be read via read.table()
#'                           otherwise, if matrix or Matrix, will use the data directly.
#' 
#' @param gene_order_file data file containing the positions of each gene along each chromosome in the genome.
#'
#' @param annotations_file a description of the cells, indicating the cell type classifications
#'
#' @param ref_group_names a vector containing the classifications of the reference (normal) cells to use for infering cnv
#'
#' @param delim delimiter used in the input files
#'
#' @param max_cells_per_group maximun number of cells to use per group. Default=NULL, using all cells defined in the annotations_file. This option is useful for randomly subsetting the existing data for a quicker preview run, such as using 50 cells per group instead of hundreds.
#' 
#' @param min_max_counts_per_cell minimum and maximum counts allowed per cell. Any cells outside this range will be removed from the counts matrix. default=(100, +Inf) and uses all cells. If used, should be set as c(min_counts, max_counts)
#'
#' @param chr_exclude list of chromosomes in the reference genome annotations that should be excluded from analysis.  Default = c('chrX', 'chrY', 'chrM')
#'
#' @description Creation of an infercnv object. This requires the following inputs:
#' A more detailed description of each input is provided below:
#'
#' The raw_counts_matrix:
#'
#'           MGH54_P16_F12 MGH53_P5_C12 MGH54_P12_C10 MGH54_P16_F02 MGH54_P11_C11  ...
#' DDX11L1     0.0000000     0.000000      0.000000      0.000000     0.0000000
#' WASH7P      0.0000000     2.231939      7.186235      5.284944     0.9650009
#' FAM138A     0.1709991     0.000000      0.000000      0.000000     0.0000000
#' OR4F5       0.0000000     0.000000      0.000000      0.000000     0.0000000
#' OR4F29      0.0000000     0.000000      0.000000      0.000000     0.0000000
#' ...
#'
#' The gene_order_file, contains chromosome, start, and stop position for each gene, tab-delimited:
#'
#'          chr  start   stop
#' DDX11L1 chr1  11869  14412
#' WASH7P  chr1  14363  29806
#' FAM138A chr1  34554  36081
#' OR4F5   chr1  69091  70008
#' OR4F29  chr1 367640 368634
#' OR4F16  chr1 621059 622053
#' ...
#' 
#' The annotations_file, containing the cell name and the cell type classification, tab-delimited.
#'
#'             V1                   V2
#' 1 MGH54_P2_C12 Microglia/Macrophage
#' 2 MGH36_P6_F03 Microglia/Macrophage
#' 3 MGH53_P4_H08 Microglia/Macrophage
#' 4 MGH53_P2_E09 Microglia/Macrophage
#' 5 MGH36_P5_E12 Oligodendrocytes (non-malignant)
#' 6 MGH54_P2_H07 Oligodendrocytes (non-malignant)
#' ...
#' 179  93_P9_H03 malignant
#' 180 93_P10_D04 malignant
#' 181  93_P8_G09 malignant
#' 182 93_P10_B10 malignant
#' 183  93_P9_C07 malignant
#' 184  93_P8_A12 malignant
#' ...
#'
#'
#' and the ref_group_names vector might look like so:  c("Microglia/Macrophage","Oligodendrocytes (non-malignant)")
#'
#' @return infercnv
#'
#' @export
#'
#' @examples
#' data(infercnv_data_example)
#' data(infercnv_annots_example)
#' data(infercnv_genes_example)
#'
#' infercnv_object_example <- infercnv::CreateInfercnvObject(raw_counts_matrix=infercnv_data_example, 
#'                                                gene_order_file=infercnv_genes_example,
#'                                                annotations_file=infercnv_annots_example,
#'                                                ref_group_names=c("normal"))
#'










CreateInfercnvObject <- function(raw_counts_matrix,
                                 gene_order_file,
                                 annotations_file,
                                 ref_group_names,
                                 delim="\t",
                                 max_cells_per_group=NULL,
                                 min_max_counts_per_cell=c(100, +Inf), # can be c(low,high) for colsums
                                 chr_exclude=c('chrX', 'chrY', 'chrM') ) {
    
    ## input expression data
    if (Reduce("|", is(raw_counts_matrix) == "character")) {
        flog.info(sprintf("Parsing matrix: %s", raw_counts_matrix)) 

        if (substr(raw_counts_matrix, nchar(raw_counts_matrix)-2, nchar(raw_counts_matrix)) == ".gz") {
            raw.data <- read.table(connection <- gzfile(raw_counts_matrix, 'rt'), sep=delim, header=TRUE, row.names=1, check.names=FALSE)
            close(connection)
            raw.data <- as.matrix(raw.data)
        }
        else if(substr(raw_counts_matrix, nchar(raw_counts_matrix)-3, nchar(raw_counts_matrix)) == ".rds") {
            raw.data <- readRDS(raw_counts_matrix)
        }
        else {
            raw.data <- read.table(raw_counts_matrix, sep=delim, header=TRUE, row.names=1, check.names=FALSE)
            raw.data <- as.matrix(raw.data)
        }
    } else if (Reduce("|", is(raw_counts_matrix) %in% c("dgCMatrix", "matrix"))) {
        # use as is:
        raw.data <- raw_counts_matrix
    } else if (Reduce("|", is(raw_counts_matrix) %in% c("data.frame"))) {
        raw.data <- as.matrix(raw_counts_matrix)
    } else {
        stop("CreateInfercnvObject:: Error, raw_counts_matrix isn't recognized as a matrix, data.frame, or filename")
    }

    ## get gene order info
    if (Reduce("|", is(gene_order_file) == "character")) {
        flog.info(sprintf("Parsing gene order file: %s", gene_order_file))
        gene_order <- read.table(gene_order_file, header=FALSE, row.names=1, sep="\t", check.names=FALSE)
    }
    else if (Reduce("|", is(gene_order_file) %in% c("dgCMatrix", "matrix", "data.frame"))) {
        gene_order <- gene_order_file
    }
    else {
        stop("CreateInfercnvObject:: Error, gene_order_file isn't recognized as a matrix, data.frame, or filename")
    }
    names(gene_order) <- c(C_CHR, C_START, C_STOP)
    if (! is.null(chr_exclude) && any(which(gene_order$chr %in% chr_exclude))) {
        gene_order = gene_order[-which(gene_order$chr %in% chr_exclude),]
    }
    
    ## read annotations file
    if (Reduce("|", is(annotations_file) == "character")) {
        flog.info(sprintf("Parsing cell annotations file: %s", annotations_file))
        input_classifications <- read.table(annotations_file, header=FALSE, row.names=1, sep=delim, stringsAsFactors=FALSE, colClasses = c('character', 'character'))
    }
    else if (Reduce("|", is(annotations_file) %in% c("dgCMatrix", "matrix", "data.frame"))) {
        input_classifications <- annotations_file
    }
    else {
        stop("CreateInfercnvObject:: Error, annotations_file isn't recognized as a matrix, data.frame, or filename")
    }



    ## just in case the first line is a default header, remove it:
    if (rownames(input_classifications)[1] == "V1") {
        input_classifications = input_classifications[-1, , drop=FALSE]
    }
    
    ## make sure all reference samples are accounted for:
    if (! all( rownames(input_classifications) %in% colnames(raw.data)) ) {
        
        missing_cells <- rownames(input_classifications)[ ! ( rownames(input_classifications) %in% colnames(raw.data) ) ]

        error_message <- paste("Please make sure that all the annotated cell ",
                               "names match a sample in your data matrix. ",
                               "Attention to: ",
                               paste(missing_cells, collapse=","))
        stop(error_message)
    }

    ## extract the genes indicated in the gene ordering file:
    order_ret <- .order_reduce(data=raw.data, genomic_position=gene_order)

    num_genes_removed = dim(raw.data)[1] - dim(order_ret$exp)[1]

    if (num_genes_removed > 0) {
        flog.info(paste("num genes removed taking into account provided gene ordering list: ",
                        num_genes_removed,
                        " = ",
                        num_genes_removed / dim(raw.data)[1] * 100,
                        "% removed.", sep=""))
    }
    
    raw.data <- order_ret$expr
    input_gene_order <- order_ret$order
    input_gene_order[["chr"]] = droplevels(input_gene_order[["chr"]])
    
    if(is.null(raw.data)) {
        error_message <- paste("None of the genes in the expression data",
                               "matched the genes in the reference genomic",
                               "position file. Analysis Stopped.")
        stop(error_message)
    }

    ## Determine if we need to do filtering on counts per cell
    if (is.null(min_max_counts_per_cell)) {
        min_max_counts_per_cell = c(1, +Inf)
    }

    min_counts_per_cell = max(1, min_max_counts_per_cell[1])  # so that it is always at least 1
    max_counts_per_cell = min_max_counts_per_cell[2]

    cs = colSums(raw.data)

    cells.keep <- which(cs >= min_counts_per_cell & cs <= max_counts_per_cell)

    n_orig_cells <- ncol(raw.data)
    n_to_remove <- n_orig_cells - length(cells.keep)
    
    flog.info(sprintf("-filtering out cells < %g or > %g, removing %g %% of cells",
                      min_counts_per_cell,
                      max_counts_per_cell,
                      n_to_remove/n_orig_cells * 100) )
    
    raw.data <- raw.data[, cells.keep]
    
    input_classifications <- input_classifications[ rownames(input_classifications) %in% colnames(raw.data), , drop=FALSE]

    orig_ref_group_names = ref_group_names
    ref_group_names <- ref_group_names[ ref_group_names %in% unique(input_classifications[,1]) ]
    if (! all.equal(ref_group_names, orig_ref_group_names)) {
        flog.warn(sprintf("-warning, at least one reference group has been removed due to cells lacking: %s",
                          orig_ref_group_names[! orig_ref_group_names %in% ref_group_names ] ))
    }

    
    
    if (! is.null(max_cells_per_group)) {
        ## trim down where needed.
        grps = split(input_classifications, input_classifications[,1])
        newdf = NULL
        for (grp in names(grps)) {
            df = grps[[grp]]
            if (dim(df)[1] > max_cells_per_group) {
                flog.info(sprintf("-reducing number of cells for grp %s from %g to %g",
                                  grp, dim(df)[1], max_cells_per_group))

                grps[[grp]] = df[sample(seq_len(dim(df)[1]), max_cells_per_group),,drop=FALSE]
            }
        }
        input_classifications = data.frame(Reduce(rbind, grps))
    }
        
    ## restrict expression data to the annotated cells.
    raw.data <- raw.data[,colnames(raw.data) %in% rownames(input_classifications)]
        
    ## reorder cell classifications according to expression matrix column names
    input_classifications <- input_classifications[order(match(row.names(input_classifications), colnames(raw.data))), , drop=FALSE]

    
    ## get indices for reference cells
    ref_group_cell_indices = list()
    for (name_group in ref_group_names) {
        cell_indices = which(input_classifications[,1] == name_group)
        
        if (length(cell_indices) == 0 ) {
            stop(sprintf("Error, not identifying cells with classification %s", name_group))
        }
        ref_group_cell_indices[[ name_group ]] <- cell_indices
    }
    
    ## rest of the cells are the 'observed' set.
    all_group_names <- unique(input_classifications[,1])
    obs_group_names <- sort(setdiff(all_group_names, ref_group_names))

    ## define groupings according to the observation annotation names
    
    obs_group_cell_indices = list()
    for (name_group in obs_group_names) {
        cell_indices = which(input_classifications[,1] == name_group)
        obs_group_cell_indices[[ toString(name_group) ]] <- cell_indices
    }

    if ((2*ncol(raw.data)) >=  10^getOption("scipen")) {
        flog.warn(paste0("Please use \"options(scipen = 100)\" before running infercnv ",
                         "if you are using the analysis_mode=\"subclusters\" option or ",
                         "you may encounter an error while the hclust is being generated."))
    }
    
    object <- new(
        Class = "infercnv",
        expr.data = raw.data, 
        count.data = raw.data,
        gene_order = input_gene_order,
        reference_grouped_cell_indices = ref_group_cell_indices,
        observation_grouped_cell_indices = obs_group_cell_indices,
        tumor_subclusters = NULL,
        options = list("chr_exclude" = chr_exclude,
                       "max_cells_per_group" = max_cells_per_group,
                       "min_max_counts_per_cell" = min_max_counts_per_cell,
                       "counts_md5" = digest(raw.data)),
        .hspike = NULL)

    validate_infercnv_obj(object)
    
    return(object)
}





























