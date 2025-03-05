#' Merging Visium spatial transciptomics count and annotation data, as well as applying a QC filter to only include spots with >= 500 counts
#'
#' MergingCountAndAnnotationData()
#' 
#' @param SectionName A character string for section name.
#' @param InputAnnotationFile An annotation file containing all barcodes to be used in the analysis (bound dataframe of one or more outputs from ImportHistologicalAnnotations())
#' @param InputCountFile A dataframe of Visium count data (output from ImportCountData())
#' @return A dataframe of barcodes with appended section names that have passed QC
#' @examples
#' MergingCountAndAnnotationData("H2_1",MergedAll, H2_1_ENSBMLID_Counts)

MergingCountAndAnnotationData <- function(SectionName, InputAnnotationFile, InputCountFile) {
  formerge <- select(InputAnnotationFile, -Histology)
  MergedAnnotationsandCounts <- inner_join(formerge, InputCountFile)
  MergedAnnotationsandCounts <- remove_rownames(MergedAnnotationsandCounts)
  MergedAnnotationsandCounts <- column_to_rownames(MergedAnnotationsandCounts, "Barcode")
  #str(MergedAnnotationsandCounts)  # 查看数据类型
  #head(MergedAnnotationsandCounts) # 查看前几行数据
  print("ok1")
  print(str(MergedAnnotationsandCounts))  # 查看数据类型
  MergedAnnotationsandCounts$Total <- rowSums(
  MergedAnnotationsandCounts[, sapply(MergedAnnotationsandCounts, is.numeric), drop = FALSE],
  na.rm = TRUE
)
  
  #MergedAnnotationsandCounts$Total <- rowSums(MergedAnnotationsandCounts)
  summary(MergedAnnotationsandCounts$Total)
  if (!"Total" %in% colnames(MergedAnnotationsandCounts)) {
  stop("Error: 'Total' column not found after rowSums() calculation!")
}
  print("ok2")
  MergedAnnotationsandCounts <- MergedAnnotationsandCounts %>% filter(MergedAnnotationsandCounts$Total >= 500)
   print("ok3")
  MergedAnnotationsandCounts <- select(MergedAnnotationsandCounts, -Total)
  MergedAnnotationsandCounts <- as.data.frame(t(MergedAnnotationsandCounts))
  MergedAnnotationsandCounts <- MergedAnnotationsandCounts[,colSums(is.na(MergedAnnotationsandCounts))<nrow(MergedAnnotationsandCounts)]
  MergedAnnotationsandCounts <- tibble::rownames_to_column(MergedAnnotationsandCounts, "Genes")
  return(MergedAnnotationsandCounts)
}
