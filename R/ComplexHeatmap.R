# Work in progress
# Function is intended to create a heatmap of the inferCNV output using the ComplexHeatmap package

# These have not yet been extensively tested and are subject to change

## To be fixed
# 1. Not plot ref plot when no ref
# 2. Clone annotation ordered with individual dendrograms
# 3. Add option for custom colors
# 4. Histology option (where does Wencheng store it?)

plot_complex_heatmap <- function(infercnv_obj,
                                 output_dir,
                                 heatmap_thresholds_file,
                                 observation_groupings_file,
                                 dendrogram_file,
                                 output_pdf_name = "complex_heatmap.pdf") {
  
  library(ComplexHeatmap)
  library(tidyverse)
  library(RColorBrewer)
  library(circlize)
  library(grid)
  library(gridExtra)
  library(ape)
  library(ggplot2)
  
  # Fetch data
  expr <- infercnv_obj@expr.data
  ref_loc <- infercnv_obj@reference_grouped_cell_indices
  obs_loc <- infercnv_obj@observation_grouped_cell_indices
  
  # Sample colors
  unique_samples <- unique(sub("^(.*?)_.*", "\\1", colnames(expr)))
  if (length(unique_samples) <= 9) {
    sample_colors <- setNames(RColorBrewer::brewer.pal(length(unique_samples), "Set1"), unique_samples)
  } else {
    sample_colors <- setNames(
      colorRampPalette(RColorBrewer::brewer.pal(9, "Set1"))(length(unique_samples)), unique_samples)
  }
  sample_colors <- sample_colors[!is.na(names(sample_colors))]
  
  # Gene ordering
  gn <- rownames(expr)
  geneOrder <- infCNV_out@gene_order
  geneOrder <- geneOrder %>% filter(!str_detect(chr, "chrM"))
  sub_geneOrder <- geneOrder[intersect(gn, rownames(geneOrder)), ]
  expr <- expr[intersect(gn, rownames(geneOrder)), ]

  # Chromosome colors
  chr_colors <- colorRampPalette(RColorBrewer::brewer.pal(12, "Set3"))(length(unique(geneOrder$chr)))
  
  # CNV Color function
  cnv_colors <- c("#00008B","#24249B","#4848AB","#6D6DBC","#9191CC","#B6B6DD","#DADAEE","#FFFFFF",
                  "#EEDADA","#DDB6B6","#CC9191","#BC6D6D","#AB4848","#9B2424","#8B0000")
  cnv_breakpoints <- read_csv(heatmap_thresholds_file, col_names = FALSE)[[1]]
  cnv_col_fun <- function(x) {
    idx <- findInterval(x, cnv_breakpoints)
    idx[idx < 1] <- 1
    idx[idx > length(cnv_colors)] <- length(cnv_colors)
    cnv_colors[idx]
  }
  
  # Histology colors
  hist_colors_default <- c("#0077B6","#F4C0BA","#F7952B","#E54C2E","#DD0303","#7C2323","#FFEA00","#29BF12","#212529","#6C757D","#AFE8F4")
  
  # Class colors
  class_colors_default <- c("#CCCCCC","#736F72","#CAF0F8","#80FFDB","#00B4D8","#01497C","#155D27","#29BF12","#ABFF4F",
                            "#FCF300","#FF8000","#E54C2E","#D40000","#780000","#9B3F11","#A24CCD","#5A189A","#000000")
  
  # Annotations - Reference
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
  
  ref_order <- infercnv_obj@tumor_subclusters$hc[[2]]$labels[infercnv_obj@tumor_subclusters$hc[[2]]$order]
  ref_expr <- ref_expr[, ref_order]
  ref_anno_df <- ref_anno_df[ref_order, ]
  
  # Annotations - Observations
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
  
  # Observation grouping
  obs_grouping <- read.table(observation_groupings_file, header = TRUE, row.names = 1, sep = " ", stringsAsFactors = FALSE)
  obs_grouping$class <- obs_anno_df[rownames(obs_grouping), "class"]
  
  ###################### Just for testing, remove once histology is added
  obs_anno_df$histology <- obs_anno_df$class
  obs_anno_df$class <- c(
    rep("Clone_A", 390),
    rep("Clone_B", 390),
    rep("Clone_C", 390),
    rep("Clone_D", 1563 - 3 * 390)
  )
  ########################################################################
  
  # Obs annotation colors
  hist_colors <- setNames(hist_colors_default, unique(obs_anno_df$histology))
  hist_colors <- hist_colors[!is.na((names(hist_colors)))]
  class_levels <- sort(unique(obs_anno_df$class))
  class_color_indices <- round(seq(1, length(class_colors_default), length.out = length(class_levels)))
  class_colors <- setNames(class_colors_default[class_color_indices], class_levels)
  
  # Dendrogram
  tree_list <- read.tree(dendrogram_file)
  dend <- as.dendrogram(as.hclust(tree_list))

  obs_expr <- obs_expr[, labels(dend)]
  obs_anno_df <- obs_anno_df[labels(dend), ]
  
  # Set ht_options
  ht_opt$ROW_ANNO_PADDING = unit(3, "mm")
  ht_opt$DENDROGRAM_PADDING = unit(3, "mm")
  ht_opt$HEATMAP_LEGEND_PADDING = unit(5, "mm")
  ht_opt$ANNOTATION_LEGEND_PADDING = unit(10, "mm")
  ht_opt$message = FALSE
  
  # Heatmap objects
  ref_ht <- Heatmap(
    t(ref_expr),
    col = cnv_col_fun,
    cluster_rows = TRUE,
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
  
  top_anno <- HeatmapAnnotation(foo = anno_block(gp = gpar(fill = chr_colors, col = NA), labels = NULL))

  bottom_anno <- HeatmapAnnotation(
    foo = anno_block(gp = gpar(fill = "NA", col = "NA"),
                     labels = unique(geneOrder$chr),
                     labels_rot = 90,
                     labels_gp = gpar(cex = 1))
  )
  
  obs_left_anno <- rowAnnotation(
    df = obs_anno_df[, c("sample", "histology", "class")],
    col = list(sample = sample_colors,
    histology = hist_colors,
    class = class_colors),
    width = unit(19, "mm"),
    simple_anno_size_adjust = TRUE,
    gap = unit(3, "mm"),
    annotation_legend_param = list(
      sample = list(
        direction = "horizontal",
        nrow = 1,
        grid_height = unit(8, "mm"),
        grid_width  = unit(8, "mm"),
        labels_gp   = gpar(fontsize = 16),
        title_gp    = gpar(fontsize = 16,
        fontface = "bold")),
      class = list(
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
        fontface = "bold"))
    )
  )
  
  obs_ht <- Heatmap(
    t(obs_expr),
    col = cnv_col_fun,
    cluster_rows = rev(dend),
    cluster_columns = FALSE,
    show_column_names = FALSE,
    show_row_names = FALSE,
    column_split = factor(sub_geneOrder$chr, levels = unique(sub_geneOrder$chr)),
    column_gap = unit(0, "mm"),
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
  
  ht_list <- ref_ht %v% obs_ht
  
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
  
  # Output
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