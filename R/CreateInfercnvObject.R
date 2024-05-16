
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





























