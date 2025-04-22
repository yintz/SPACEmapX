# Work in progress
# These have not yet been extensively tested and are subject to change

#' Plot Complex Heatmap for inferCNV Output
#'
#' This function generates a heatmap using the ComplexHeatmap package to visualize inferCNV results,
#' supporting reference and observation groups, histology and clone annotations, and user-defined color schemes.
#' 
#' @param infercnv_obj An inferCNV object created by the `infercnv::create_infercnv_object` and processed by `run()` or `run_infercnv()`.
#' @param output_dir Directory where the output heatmap PDF will be saved.
#' @param heatmap_thresholds_file Path to the inferCNV-generated `heatmap_thresholds.txt` file.
#' @param dendrogram_file A Newick file for dendrogram, can be either a single tree (`phylo`) or multiple trees (`multiPhylo`).
#' @param obs_groups Character string indicating how observations are grouped: `"histology"` (default) or `"clone"`.
#' @param hist_annotation A data frame or path to a TSV file with histology information (must include a `"Histology"` column).
#' @param output_pdf_name Filename for the saved heatmap PDF. Default is `"complex_heatmap.pdf"`.
#' @param hist_cols Optional. Named list or vector of colors for histology groups. Format: `c("group1" = "#hex", ...)`.
#' @param clone_cols Optional. Named list or vector of colors for clone groups. Format: `c("clone1" = "#hex", ...)`.
#'
#' @return A PDF file is saved to `output_dir`, visualizing CNV patterns across references and observations with dendrograms and annotations.
#'
#' @details
#' - For `dendrogram_file`, a single tree (`phylo`) will apply a global dendrogram, while a multiple tree (`multiPhylo`) input will order spots per clone.
#' - The function will auto-detect and adjust for reference presence.
#' - Colors for sample, histology, and clone groups are automatically generated but can be overridden.
#' - Chromosome and CNV coloring is derived from inferCNV settings and thresholds.
#'
#' @examples
#' plot_complex_heatmap(
#'   infercnv_obj = my_infercnv,
#'   output_dir = "results/",
#'   heatmap_thresholds_file = "results/infercnv.heatmap_thresholds.txt",
#'   dendrogram_file = "results/infercnv.observations_dendrogram.txt",
#'   obs_groups = "clone",
#'   hist_annotation = "metadata/histologies.tsv",
#'   hist_cols = c("Tumor" = "#FF0000", "Normal" = "#0000FF"),
#'   clone_cols = c("Clone1" = "#E41A1C", "Clone2" = "#377EB8")
#' )

plot_complex_heatmap <- function(infercnv_obj,
                                 output_dir,
                                 heatmap_thresholds_file,
                                 dendrogram_file,
                                 obs_groups = "histology",
                                 hist_annotation,
                                 output_pdf_name = "complex_heatmap.pdf",
                                 hist_cols = NULL,
                                 clone_cols = NULL) {
  
  library(ComplexHeatmap)
  library(tidyverse)
  library(RColorBrewer)
  library(circlize)
  library(grid)
  library(gridExtra)
  library(ape)
  library(ggplot2)
  library(dendextend)
  
  if (!obs_groups %in% c("histology", "clone")) {
    stop("`obs_groups` must be either 'histology' or 'clone'.")
  }
  
  # Fetch expressions data
  expr <- infercnv_obj@expr.data
  
  # Gene ordering
  gn <- rownames(expr)
  geneOrder <- infercnv_obj@gene_order
  geneOrder <- geneOrder %>% filter(!str_detect(chr, "chrM"))
  sub_geneOrder <- geneOrder[intersect(gn, rownames(geneOrder)), ]
  expr <- expr[intersect(gn, rownames(geneOrder)), ]
  
  # Colors
  ## Sample
  unique_samples <- unique(sub("^(.*?)_.*", "\\1", colnames(expr)))
  if (length(unique_samples) <= 9) {
    sample_colors <- setNames(RColorBrewer::brewer.pal(length(unique_samples), "Set1"), unique_samples)
  } else {
    sample_colors <- setNames(
      colorRampPalette(RColorBrewer::brewer.pal(9, "Set1"))(length(unique_samples)), unique_samples)
  }
  sample_colors <- sample_colors[!is.na(names(sample_colors))]
  
  ## Chromosome
  chr_colors <- colorRampPalette(RColorBrewer::brewer.pal(12, "Set3"))(length(unique(geneOrder$chr)))
  
  ## CNV
  cnv_colors <- c("#00008B","#24249B","#4848AB","#6D6DBC","#9191CC","#B6B6DD","#DADAEE","#FFFFFF",
                  "#EEDADA","#DDB6B6","#CC9191","#BC6D6D","#AB4848","#9B2424","#8B0000")
  
  ## Histology
  hist_colors_default <- c("#0077B6","#F4C0BA","#F7952B","#E54C2E","#DD0303","#7C2323","#FFEA00","#29BF12","#212529","#6C757D","#AFE8F4")
  
  ## Clone
  clone_colors_default <- c("#CCCCCC","#736F72","#CAF0F8","#80FFDB","#00B4D8","#01497C","#155D27","#29BF12","#ABFF4F",
                            "#FCF300","#FF8000","#E54C2E","#D40000","#780000","#9B3F11","#A24CCD","#5A189A","#000000")
  
  cnv_breakpoints <- read_csv(heatmap_thresholds_file, col_names = FALSE, show_col_types = FALSE)[[1]]
  cnv_col_fun <- function(x) {
    idx <- findInterval(x, cnv_breakpoints)
    idx[idx < 1] <- 1
    idx[idx > length(cnv_colors)] <- length(cnv_colors)
    cnv_colors[idx]
  }
  
  # Annotations - Observations
  obs_loc <- infercnv_obj@observation_grouped_cell_indices
  names(obs_loc) <- make.names(names(obs_loc))
  list2env(obs_loc, envir = environment())
  
  obs_anno_list <- list()
  for (annotation in names(obs_loc)) {
    obs_anno_list[[annotation]] <- data.frame(
      CB = colnames(expr)[get(annotation)],
      class = annotation
    )
  }
  obs_anno_df <- do.call(rbind, obs_anno_list)
  
  obs_expr <- expr[, colnames(expr) %in% obs_anno_df$CB]
  obs_anno_df <- obs_anno_df[match(colnames(obs_expr), obs_anno_df$CB), ]
  
  rownames(obs_anno_df) <- obs_anno_df$CB
  obs_anno_df$CB <- NULL
  
  obs_anno_df$sample <- sub("^(.*?)_.*", "\\1", rownames(obs_anno_df))
  
  # Dendrogram
  tree <- read.tree(dendrogram_file)
  if (inherits(tree, "phylo")) {
    dend <- as.dendrogram(as.hclust(tree))
    dend <- rev(dend)
    
    # Re-order based on dendrogram
    dend_labels <- labels(dend)
    missing_labels <- setdiff(dend_labels, colnames(obs_expr))
    if (length(missing_labels) > 0) {
      stop("Some labels in the dendrogram are not found in obs_expr: ", paste(missing_labels, collapse = ", "))
    }
    obs_expr <- obs_expr[, labels(dend)]
    obs_anno_df <- obs_anno_df[labels(dend), ]
    
    row_split <- NULL
  }
  else if (inherits(tree, "multiPhylo")) {
    
    clone_tree_list <- read.tree(dendrogram_file)
    
    clone_dend_list <- lapply(clone_tree_list, function(tree) {
      as.dendrogram(as.hclust(tree))
    })
    
    clone_order_list <- lapply(clone_dend_list, labels)
    clone_spot_order <- unlist(clone_order_list)
    names(clone_order_list) <- names(infercnv_obj@observation_grouped_cell_indices)
    
    clone_names <- names(clone_order_list)
    clone_lengths <- sapply(clone_order_list, length)
    row_split <- factor(rep(clone_names, times = clone_lengths), levels = clone_names) 
    
    # Re-order based on dendrogram
    missing_labels <- setdiff(clone_spot_order, colnames(obs_expr))
    if (length(missing_labels) > 0) {
      stop("Some labels in the dendrogram are not found in obs_expr: ", paste(missing_labels, collapse = ", "))
    }
    obs_expr <- obs_expr[, clone_spot_order]
    obs_anno_df <- obs_anno_df[clone_spot_order, ]
    
    dend <- function(mat) hclust(dist(mat), method = "ward.D2")
  }
  
  # Set ht_options
  ht_opt$ROW_ANNO_PADDING = unit(3, "mm")
  ht_opt$DENDROGRAM_PADDING = unit(3, "mm")
  ht_opt$HEATMAP_LEGEND_PADDING = unit(5, "mm")
  ht_opt$ANNOTATION_LEGEND_PADDING = unit(10, "mm")
  ht_opt$message = FALSE
  
  # Obs Heatmap object
  top_anno <- HeatmapAnnotation(foo = anno_block(gp = gpar(fill = chr_colors, col = NA), 
                                                 labels = NULL), 
                                show_annotation_name = FALSE)
  
  bottom_anno <- HeatmapAnnotation(
    foo = anno_block(gp = gpar(fill = "NA", 
                               col = "NA"),
                     labels = unique(geneOrder$chr),
                     labels_rot = 90,
                     labels_gp = gpar(cex = 1))
  )
  
  if (obs_groups == "histology"){
    obs_anno_df <- obs_anno_df %>% rename(histology = class)
    
    if (!is.null(hist_cols)) {
      hist_colors <- hist_cols
    } else {
      hist_colors <- setNames(hist_colors_default, 
                              sort(unique(obs_anno_df$histology)))
    }
    
    hist_colors <- hist_colors[!is.na((names(hist_colors)))]
    
    anno_type <- c("sample", "histology")
    anno_col <- list(sample = sample_colors,
                     histology = hist_colors)
    anno_leg <- list(
      sample = list(
        direction = "horizontal",
        nrow = 1,
        grid_height = unit(8, "mm"),
        grid_width  = unit(8, "mm"),
        labels_gp   = gpar(fontsize = 16),
        title_gp    = gpar(fontsize = 16,
                           fontface = "bold")
      ),
      histology = list(
        direction = "horizontal",
        nrow = 1,
        grid_height = unit(8, "mm"),
        grid_width  = unit(8, "mm"),
        labels_gp   = gpar(fontsize = 16),
        title_gp    = gpar(fontsize = 16,
                           fontface = "bold")
      )
    )
    
    obs_left_anno <- rowAnnotation(
      df = obs_anno_df[, anno_type],
      col = anno_col,
      width = unit(19, "mm"),
      simple_anno_size_adjust = TRUE,
      gap = unit(3, "mm"),
      annotation_legend_param = anno_leg
    )
  }
  else if (obs_groups == "clone"){
    obs_anno_df <- obs_anno_df %>% rename(clone = class)
    
    if (!is.null(clone_cols)) {
      clone_colors <- clone_cols
    } else {
      clone_levels <- sort(unique(obs_anno_df$clone))
      clone_color_indices <- round(seq(1, length(clone_colors_default), length.out = length(clone_levels)))
      clone_colors <- setNames(clone_colors_default[clone_color_indices], clone_levels)
    }
    
    if (is.character(hist_annotation) && file.exists(hist_annotation)) {
      hist_anno <- read.table(hist_annotation, row.names = 1,sep = "\t", stringsAsFactors = FALSE)
      
      if (!"Histology" %in% colnames(hist_anno)) {
        stop("hist_annotation file must contain a column named 'Histology'")
      }
      
      obs_anno_df$histology <- hist_anno[rownames(obs_anno_df), "Histology"]
      
    } else if (is.data.frame(hist_annotation)) {
      
      if (!all(rownames(obs_anno_df) %in% rownames(hist_annotation))) {
        stop("hist_annotation data frame rownames must match obs_anno_df rownames")
      }
      
      obs_anno_df$histology <- hist_annotation[rownames(obs_anno_df), 1]
      
    } else {
      
      stop("`hist_annotation` must be either a file path to a .tsv or a data frame.")
      
    }
    
    if (!is.null(hist_cols)) {
      hist_colors <- hist_cols
    } else {
      hist_colors <- setNames(hist_colors_default, 
                              sort(unique(obs_anno_df$histology)))
    }
    
    hist_colors <- hist_colors[!is.na((names(hist_colors)))]
    
    anno_type <- c("sample", "histology", "clone")
    anno_col <- list(sample = sample_colors,
                     histology = hist_colors,
                     clone = clone_colors)
    anno_leg <- list(
      sample = list(
        direction = "horizontal",
        nrow = 1,
        grid_height = unit(8, "mm"),
        grid_width  = unit(8, "mm"),
        labels_gp   = gpar(fontsize = 16),
        title_gp    = gpar(fontsize = 16,
                           fontface = "bold")),
      clone = list(
        direction = "horizontal",
        nrow = 1,
        grid_height = unit(8, "mm"),
        grid_width  = unit(8, "mm"),
        labels_gp   = gpar(fontsize = 16),
        title_gp    = gpar(fontsize = 16,
                           fontface = "bold")),
      histology = list(
        direction = "horizontal",
        nrow = 1,
        grid_height = unit(8, "mm"),
        grid_width  = unit(8, "mm"),
        labels_gp   = gpar(fontsize = 16),
        title_gp    = gpar(fontsize = 16,
                           fontface = "bold")
      )
    )
    
    obs_left_anno <- rowAnnotation(
      df = obs_anno_df[, anno_type],
      col = anno_col,
      width = unit(30, "mm"),
      simple_anno_size_adjust = TRUE,
      gap = unit(3, "mm"),
      annotation_legend_param = anno_leg
    )
  }
  
  obs_ht <- Heatmap(
    t(obs_expr),
    col = cnv_col_fun,
    cluster_rows = dend,
    cluster_columns = FALSE,
    show_column_names = FALSE,
    show_row_names = FALSE, 
    show_row_dend = TRUE,
    row_dend_side = "left",
    column_split = factor(sub_geneOrder$chr, levels = unique(sub_geneOrder$chr)),
    column_gap = unit(0, "mm"), 
    row_split = row_split, 
    cluster_row_slices = FALSE, 
    row_gap = unit(0, "mm"),
    border = TRUE,
    show_heatmap_legend = FALSE,
    bottom_annotation = bottom_anno,
    top_annotation = top_anno,
    left_annotation = obs_left_anno,
    row_title = "Observations (Spots)",
    row_title_side = "right",
    row_title_rot = 90,
    row_title_gp = gpar(fontsize = 18, fontface = "bold", lineheight = 3),
    row_dend_width = unit(6, "cm"),
    height = unit(8, "inch")
  )
  
  # Function for heatmap legend
  plot_infercnv_legend_static <- function(expr_data,
                                          cnv_breakpoints,
                                          cnv_colors,
                                          denscol = "blue",
                                          densadj = 0.25,
                                          key.xlab = "Modified Expression") {
    
    round_limits <- function(x) {
      values <- sort(unique(c(seq(0, 10, 0.10), seq(0, 10, 0.15))))
      nearest <- values[which.min(abs(values - x))]
      return(nearest)
    }
    
    x.range <- c(round_limits(cnv_breakpoints[1]), round_limits(cnv_breakpoints[length(cnv_breakpoints)]))
    
    breaks <- cnv_breakpoints
    colors <- cnv_colors
    
    expr_data_clamped <- pmax(expr_data, min(breaks))
    expr_data_clamped <- pmin(expr_data_clamped, max(breaks))
    
    h <- hist(expr_data_clamped, breaks = breaks, plot = FALSE)
    hist_data <- data.frame(
      x = rep(h$breaks, each = 2),
      y = c(0, rep(h$counts, each = 2), 0)
    )
    hist_data$y <- hist_data$y / max(hist_data$y, na.rm = TRUE) * 0.98
    
    color_data <- data.frame(
      xstart   = breaks[-length(breaks)],
      xend     = breaks[-1],
      midpoint = 0.5 * (breaks[-length(breaks)] + breaks[-1])
    )
    
    p <- ggplot() +
      geom_rect(data = color_data,
                aes(xmin = xstart, xmax = xend, ymin = 0, ymax = 1, fill = midpoint),
                color = NA) +
      scale_fill_gradientn(colors = colors,
                           limits = x.range,
                           breaks = c(x.range[1], mean(x.range), x.range[2]),
                           labels = round(c(x.range[1], mean(x.range), x.range[2]), 2),
                           guide = "none") +
      coord_cartesian(expand = FALSE) +
      xlab(key.xlab) + ylab("") +
      geom_step(data = hist_data, aes(x = x, y = y), direction = "hv", color = denscol) +
      ggtitle("Distribution of Expression") +
      theme_minimal(base_size = 12) +
      theme(axis.text.y = element_blank(),
            axis.ticks.y = element_blank(),
            plot.title = element_text(hjust = 0.5),
            legend.position = "none")
    
    return(p)
  }
  
  # Heatmap legend
  legend_plot <- plot_infercnv_legend_static(expr, cnv_breakpoints, cnv_colors)
  legend_plot_grob <- ggplotGrob(legend_plot)
  
  # Annotations - Reference
  ## if there is reference
  if (is_empty(infercnv_obj@reference_grouped_cell_indices) == FALSE){
    ref_loc <- infercnv_obj@reference_grouped_cell_indices
    
    names(ref_loc) <- make.names(names(ref_loc))
    list2env(ref_loc, envir = environment())
    
    ref_anno_list <- list()
    for (annotation in names(ref_loc)) {
      ref_anno_list[[annotation]] <- data.frame(
        CB = colnames(expr)[get(annotation)],
        class = annotation
      )
    }
    ref_anno_df <- do.call(rbind, ref_anno_list)
    
    ref_expr <- expr[, colnames(expr) %in% ref_anno_df$CB]
    ref_anno_df <- ref_anno_df[match(colnames(ref_expr), ref_anno_df$CB), ]
    
    rownames(ref_anno_df) <- ref_anno_df$CB
    ref_anno_df$CB <- NULL
    
    ref_anno_df$sample <- sub("^(.*?)_.*", "\\1", rownames(ref_anno_df))
    
    ref_order <- infercnv_obj@tumor_subclusters$hc[[length(infercnv_obj@tumor_subclusters$hc)]]$labels[infercnv_obj@tumor_subclusters$hc[[length(infercnv_obj@tumor_subclusters$hc)]]$order]
    ref_expr <- ref_expr[, ref_order]
    ref_anno_df <- ref_anno_df[ref_order, ]
    
    # Ref heatmap object
    ref_ht <- Heatmap(
      t(ref_expr),
      col = cnv_col_fun,
      cluster_rows = function(mat) hclust(dist(mat), method = "ward.D2"),
      show_row_dend = FALSE,
      cluster_columns = FALSE,
      show_column_names = FALSE,
      show_row_names = FALSE,
      column_split = factor(sub_geneOrder$chr, levels = unique(sub_geneOrder$chr)),
      column_gap = unit(0, "mm"),
      column_title = NULL,
      border = TRUE,
      row_title = "References (Spots)",
      row_title_side = "right",
      row_title_rot = 90,
      row_title_gp = gpar(fontsize = 18, fontface = "bold", lineheight = 3),
      show_heatmap_legend = FALSE,
      height = unit(4, "inch")
    )
  }
  else {
    ref_ht <- Heatmap(
      matrix(NA, nrow = 1, ncol = nrow(obs_expr)),
      cluster_rows = FALSE,
      cluster_columns = FALSE,
      show_heatmap_legend = FALSE,
      show_row_names = FALSE,
      show_column_names = FALSE,
      column_split = factor(sub_geneOrder$chr, levels = unique(sub_geneOrder$chr)),
      column_gap = unit(0, "mm"),
      column_title = NULL,
      border = FALSE,
      rect_gp = gpar(fill = NA, col = NA),
      height = unit(4, "inch")
    )
  }
  
  ht_list <- ref_ht %v% obs_ht
  
  # Plot
  pdf(file.path(output_dir, output_pdf_name), width = 20, height = 16)
  draw(ht_list,
       padding = unit(c(1, 1, 2, 1), "cm"),
       column_title = "Genomic Region",
       column_title_side = "bottom",
       column_title_gp = gpar(fontsize = 18, fontface = "bold", lineheight = 3),
       ht_gap = unit(1, "cm"),
       annotation_legend_side = "bottom"
  )
  vp <- viewport(x = 0.01, y = 0.94, width = 0.15, height = 0.25, just = c("left", "top"))
  pushViewport(vp)
  grid.draw(legend_plot_grob)
  popViewport()
  dev.off()
}