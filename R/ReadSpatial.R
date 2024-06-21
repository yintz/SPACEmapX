# Load necessary library
library(readr)

# Define the ReadSpatial function
ReadSpatial <- function(file_path) {
  # Read the CSV file into a data frame
  samples_data <- read_csv(file_path)
  
  # Check if required columns are present
  required_columns <- c("samplesID", "GeneFile", "annotation", "type")
  missing_columns <- setdiff(required_columns, colnames(samples_data))
  
  if (length(missing_columns) > 0) {
    stop(paste("The following required columns are missing:", paste(missing_columns, collapse = ", ")))
  }

# Check the number of samples
  num_samples <- nrow(samples_data)
  cat("Number of samples:", num_samples, "\n")
  
# Print the first two samples
  if (num_samples = 0) {
    cat("First sample:\n")
    print(there is no data to read)
  }
  
  if (num_samples > 1) {
    NumberCheck<-1
    while (NumberCheck<=num_samples){
      ENSBMLID_Counts[[NumberCheck]]<-ImportCountData(samples_data$samplesID[NumberCheck,1], samples_data$GeneFile[NumberCheck,2])
      Histology[[NumberCheck]] <- ImportHistologicalAnnotations(samples_data$samplesID[NumberCheck,1], samples_data$GeneFile[NumberCheck,3])
      TissueType[[NumberCheck]] <- read_csv(samples_data$samplesID[NumberCheck,1], samples_data$GeneFile[NumberCheck,4])
      NumberCheck<-NumberCheck+1
  }
  
  # Return the data frame
  return(samples_data)
}

# Example usage
# file_path <- "path/to/your/file.csv"
# samples_data <- ReadSpatial(file_path)
# head(samples_data)
