# Work in progress - not yet made into a function
# Name of the script file is temporary and will be changed in the future
# Function is intended to create a heatmap of the inferCNV output using the ComplexHeatmap package

# Currently a collection of scripts to be used for plotting functions
# These have not yet been extensively tested and are subject to change
# Scripts currently for reference free
# Reference heatmap to be added





library(Seurat)
library(tidyverse)
library(SpatialInferCNV)
library(ComplexHeatmap)
library(circlize)
library(gridtext)
library(ape)

####################################################################################################
# Plot using ComplexHeatmap denoised
####################################################################################################

## Fetch data
expr <- infCNV_noref_no_exclude@expr.data
obs_loc <- infCNV_noref_no_exclude@observation_grouped_cell_indices

## Create annotation dataframe
Benign <- obs_loc$Benign
Cancer <- obs_loc$Cancer
Lymph <- obs_loc$Lymphocytes
SV <- obs_loc$`Seminal vesicle`
Stroma <- obs_loc$Stroma

anno.df=data.frame(
  CB=c(colnames(expr)[Benign], 
       colnames(expr)[Cancer],
       colnames(expr)[Lymph],
       colnames(expr)[SV],
       colnames(expr)[Stroma]),
  class=c(rep("Benign",length(Benign)),
          rep("Cancer",length(Cancer)),
          rep("Lymphocytes",length(Lymph)),
          rep("Seminal vesicle",length(SV)),
          rep("Stroma",length(Stroma)))
)

## Match order of anno.df and expr matrix and
anno.df <- anno.df[match(colnames(expr), anno.df$CB), ]
## Add rownames to anno.df
rownames(anno.df) <- colnames(expr)
## Add sample info based on rownames
anno.df$sample <- sub("^(.*?)_.*", "\\1", rownames(anno.df))

## Import the grouping from inferCNV output
obs_grouping <- read.table("out_3_sections_NoExclude_rename_no_grouping/infercnv.observation_groupings.txt",
                 header = TRUE,
                 row.names = 1,
                 sep = " ",
                 stringsAsFactors = FALSE)

## Add class info to the obs_grouping to be able to retrieve the colors
obs_grouping$class <- anno.df[rownames(obs_grouping), "class"]
class_colors <- setNames(unique(obs_grouping$Annotation.Color), unique(obs_grouping$class))

## Generate sample colors
unique_samples <- unique(anno.df$sample)
if(length(unique_samples) <= 9) {
  sample_colors <- setNames(RColorBrewer::brewer.pal(n = length(unique_samples), name = "Set1"), unique_samples)
} else {
  sample_colors <- setNames(colorRampPalette(RColorBrewer::brewer.pal(9, "Set1"))(length(unique_samples)), unique_samples)
}

## Order genes
gn <- rownames(expr)
geneFile <- read.table("/srv/home/mengxiao.he/inferCNV/siCNV_GeneOrderFile.tsv",header = F,sep = "\t",stringsAsFactors = F)
rownames(geneFile)=geneFile$V1
geneFile <- geneFile %>% filter(!str_detect(V2, "chrM"))
sub_geneFile <-  geneFile[intersect(gn,geneFile$V1),]
expr=expr[intersect(gn,geneFile$V1),]

## Generate colors for the chromosomes
get_group_color_palette <- function(){
    return(colorRampPalette(RColorBrewer::brewer.pal(12,"Set3")))
}
chr_colors <- get_group_color_palette()(length(unique(geneFile$V2)))

## Heatmapp colors based on inferCNV output
cnv_colors <- c(
  "#00008B","#24249B","#4848AB","#6D6DBC","#9191CC",
  "#B6B6DD","#DADAEE","#FFFFFF","#EEDADA","#DDB6B6",
  "#CC9191","#BC6D6D","#AB4848","#9B2424","#8B0000"
)
cnv_breakpoints <- read_csv("out_3_sections_NoExclude_rename_no_grouping/infercnv.heatmap_thresholds.txt", col_names = FALSE)[[1]]

cnv_col_fun <- function(x) {
  idx <- findInterval(x, cnv_breakpoints)
  idx[idx < 1] <- 1
  idx[idx > length(cnv_colors)] <- length(cnv_colors)
  cnv_colors[idx]
}

## Import dendrogram from inferCNV output and format
tree_list <- read.tree("out_3_sections_NoExclude_rename_no_grouping/infercnv.observations_dendrogram.txt")
phylo_obj <- if (is.list(tree_list)) tree_list[[1]] else tree_list
dend <- as.dendrogram(as.hclust(tree_list))

## Reorder expression matrix, annotations
expr <- expr[, labels(dend)]
anno.df <- anno.df[labels(dend),]

top_anno <- HeatmapAnnotation(
  foo = anno_block(
    gp = gpar(fill = chr_colors, col = "black"),
    labels = NULL
  )
)

bottom_anno <- HeatmapAnnotation(
  foo = anno_block(
    gp = gpar(fill = "NA", col="NA"), 
    labels = unique(geneFile$V2), 
    labels_rot = 90,
    labels_gp = gpar(cex = 1)
  )
)

left_anno <- rowAnnotation(df = anno.df[, c("sample", "class")],
                           col = list(class = class_colors,
                                      sample = sample_colors), 
                           width = unit(20, "mm"), 
                           simple_anno_size_adjust = TRUE,
                           gap = unit(3, "mm"))

ht_opt$ROW_ANNO_PADDING = unit(3, "mm")
ht_opt$DENDROGRAM_PADDING = unit(3, "mm")
ht_opt$HEATMAP_LEGEND_PADDING = unit(5, "mm")
ht_opt$ANNOTATION_LEGEND_PADDING = unit(5, "mm")

pdf("out_3_sections_NoExclude_rename_no_grouping/complex_heatmap_denoise.pdf", width = 20, height = 10)
ht <- Heatmap(
  t(expr),
  col = cnv_col_fun,
  cluster_rows = rev(dend),
  cluster_columns = FALSE,
  show_column_names = FALSE,
  show_row_names = FALSE,
  column_split = factor(sub_geneFile$V2, 
                        levels = unique(sub_geneFile$V2)),
  column_gap = unit(0, "mm"), 
  border = TRUE,
  heatmap_legend_param = list(
    title = "Modified Expression",
    direction = "vertical",
    title_position = "leftcenter-rot",
    at = c(0.85, 1, 1.15),
    legend_height = unit(3, "cm")
  ), 
  bottom_annotation = bottom_anno,
  top_annotation = top_anno,
  left_annotation = left_anno,
  row_title = NULL,
  column_title = "Genomic Region", 
  column_title_side = "bottom", 
  row_dend_width = unit(6, "cm"),
)
draw(ht, 
     heatmap_legend_side = "right",
     padding = unit(c(1, 1, 1, 1), "cm")
)
dev.off()

####################################################################################################
# Plot using ComplexHeatmap pre-denoise
####################################################################################################

## Load prelimminary infercnv object
infCNV_noref_no_exclude_pre <- readRDS("out_3_sections_NoExclude_rename_no_grouping/preliminary.infercnv_obj")

# Fetch data
expr <- infCNV_noref_no_exclude_pre@expr.data
obs_loc <- infCNV_noref_no_exclude_pre@observation_grouped_cell_indices

# Create annotation dataframe
Benign <- obs_loc$Benign
Cancer <- obs_loc$Cancer
Lymph <- obs_loc$Lymphocytes
SV <- obs_loc$`Seminal vesicle`
Stroma <- obs_loc$Stroma

anno.df=data.frame(
  CB=c(colnames(expr)[Benign], 
       colnames(expr)[Cancer],
       colnames(expr)[Lymph],
       colnames(expr)[SV],
       colnames(expr)[Stroma]),
  class=c(rep("Benign",length(Benign)),
          rep("Cancer",length(Cancer)),
          rep("Lymphocytes",length(Lymph)),
          rep("Seminal vesicle",length(SV)),
          rep("Stroma",length(Stroma)))
)

# Match order of anno.df and expr matrix and
anno.df <- anno.df[match(colnames(expr), anno.df$CB), ]
# Add rownames to anno.df
rownames(anno.df) <- colnames(expr)
# Add sample info based on rownames
anno.df$sample <- sub("^(.*?)_.*", "\\1", rownames(anno.df))

# Import the grouping from inferCNV output
obs_grouping <- read.table("out_3_sections_NoExclude_rename_no_grouping/infercnv.preliminary.observation_groupings.txt",
                 header = TRUE,
                 row.names = 1,
                 sep = " ",
                 stringsAsFactors = FALSE)

# Add class info to the obs_grouping to be able to retrieve the colors
obs_grouping$class <- anno.df[rownames(obs_grouping), "class"]
class_colors <- setNames(unique(obs_grouping$Annotation.Color), unique(obs_grouping$class))

# Generate sample colors
unique_samples <- unique(anno.df$sample)
if(length(unique_samples) <= 9) {
  sample_colors <- setNames(RColorBrewer::brewer.pal(n = length(unique_samples), name = "Set1"), unique_samples)
} else {
  sample_colors <- setNames(colorRampPalette(RColorBrewer::brewer.pal(9, "Set1"))(length(unique_samples)), unique_samples)
}

# Order genes
gn <- rownames(expr)
geneFile <- read.table("/srv/home/mengxiao.he/inferCNV/siCNV_GeneOrderFile.tsv",header = F,sep = "\t",stringsAsFactors = F)
rownames(geneFile)=geneFile$V1
geneFile <- geneFile %>% filter(!str_detect(V2, "chrM"))
sub_geneFile <-  geneFile[intersect(gn,geneFile$V1),]
expr=expr[intersect(gn,geneFile$V1),]

# Generate colors for the chromosomes
get_group_color_palette <- function(){
    return(colorRampPalette(RColorBrewer::brewer.pal(12,"Set3")))
}
chr_colors <- get_group_color_palette()(length(unique(geneFile$V2)))

# Heatmapp colors based on inferCNV output
cnv_colors <- c(
  "#00008B","#24249B","#4848AB","#6D6DBC","#9191CC",
  "#B6B6DD","#DADAEE","#FFFFFF","#EEDADA","#DDB6B6",
  "#CC9191","#BC6D6D","#AB4848","#9B2424","#8B0000"
)
cnv_breakpoints <- read_csv("out_3_sections_NoExclude_rename_no_grouping/infercnv.preliminary.heatmap_thresholds.txt", col_names = FALSE)[[1]]

cnv_col_fun <- function(x) {
  idx <- findInterval(x, cnv_breakpoints)
  idx[idx < 1] <- 1
  idx[idx > length(cnv_colors)] <- length(cnv_colors)
  cnv_colors[idx]
}

# Import dendrogram from inferCNV output and format
tree_list <- read.tree("out_3_sections_NoExclude_rename_no_grouping/infercnv.preliminary.observations_dendrogram.txt")
phylo_obj <- if (is.list(tree_list)) tree_list[[1]] else tree_list
dend <- as.dendrogram(as.hclust(tree_list))

# Reorder expression matrix, annotations
expr <- expr[, labels(dend)]
anno.df <- anno.df[labels(dend),]

top_anno <- HeatmapAnnotation(
  foo = anno_block(
    gp = gpar(fill = chr_colors, col = "black"),
    labels = NULL
  )
)

bottom_anno <- HeatmapAnnotation(
  foo = anno_block(
    gp = gpar(fill = "NA", col="NA"), 
    labels = unique(geneFile$V2), 
    labels_rot = 90,
    labels_gp = gpar(cex = 1)
  )
)

left_anno <- rowAnnotation(df = anno.df[, c("sample", "class")],
                           col = list(class = class_colors,
                                      sample = sample_colors), 
                           width = unit(20, "mm"), 
                           simple_anno_size_adjust = TRUE,
                           gap = unit(3, "mm"))

ht_opt$ROW_ANNO_PADDING = unit(3, "mm")
ht_opt$DENDROGRAM_PADDING = unit(3, "mm")
ht_opt$HEATMAP_LEGEND_PADDING = unit(5, "mm")
ht_opt$ANNOTATION_LEGEND_PADDING = unit(5, "mm")

pdf("out_3_sections_NoExclude_rename_no_grouping/complex_heatmap_prelim.pdf", width = 20, height = 10)
ht <- Heatmap(
  t(expr),
  col = cnv_col_fun,
  cluster_rows = rev(dend),
  cluster_columns = FALSE,
  show_column_names = FALSE,
  show_row_names = FALSE,
  column_split = factor(sub_geneFile$V2, 
                        levels = unique(sub_geneFile$V2)),
  column_gap = unit(0, "mm"), 
  border = TRUE,
  heatmap_legend_param = list(
    title = "Modified Expression",
    direction = "vertical",
    title_position = "leftcenter-rot",
    at = c(0.85, 1, 1.15),
    legend_height = unit(3, "cm")
  ), 
  bottom_annotation = bottom_anno,
  top_annotation = top_anno,
  left_annotation = left_anno,
  row_title = NULL,
  column_title = "Genomic Region", 
  column_title_side = "bottom", 
  row_dend_width = unit(6, "cm"),
)
draw(ht, 
     heatmap_legend_side = "right",
     padding = unit(c(1, 1, 1, 1), "cm")
)
dev.off()