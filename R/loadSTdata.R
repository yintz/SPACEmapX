loadSTdata <- function(kk) {
  # check it has more than 3 columns
  if (ncol(kk) < 3) {
    stop("Error：dataframe needs at least 3 columns sampleID, ExpressionMatrixLocation and HistologyFileLocation）！")
  }

  # run through all row to read data.
  for (i in 1:nrow(kk)) {
    sample_id <- kk[i, 1]  # sampleID）
    matrix_path <- kk[i, 2]  # (Expression Path）
    histology_path <- kk[i, 3]  # (Pathology Path）

    # Load data into global environment
    file_name <- paste0(sample_id, "_matrix")
    file_name_prefix <- ImportCountData(sample_id, matrix_path)
    assign(file_name, file_name_prefix, envir = .GlobalEnv)

    file_hist <- paste0(sample_id, "_histology")
    file_hist_prefix <- ImportHistologicalAnnotations(sample_id, histology_path)
    assign(file_hist, file_hist_prefix, envir = .GlobalEnv)
  }
  
  print("successful！🚀")
}
