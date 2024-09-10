
# 
# # Function to sum 6x6 blocks of an 837x837 matrix
# sum_blocks <- function(mat, block_size) {
#   # Get the dimensions of the input matrix
#   nrow <- dim(mat)[1]
#   ncol <- dim(mat)[2]
#   
#   # Check if the dimensions are divisible by block_size
#   if (nrow %% block_size != 0 | ncol %% block_size != 0) {
#     stop("Matrix dimensions must be divisible by block size.")
#   }
#   
#   # Calculate the dimensions of the resulting matrix
#   new_nrow <- nrow / block_size
#   new_ncol <- ncol / block_size
#   
#   # Initialize the resulting matrix
#   result <- matrix(0, nrow = new_nrow, ncol = new_ncol)
#   
#   # Sum blocks of the matrix
#   for (i in 1:new_nrow) {
#     for (j in 1:new_ncol) {
#       # Define the indices for the current block
#       row_indices <- ((i - 1) * block_size + 1):(i * block_size)
#       col_indices <- ((j - 1) * block_size + 1):(j * block_size)
#       
#       # Sum the block and store it in the result matrix
#       result[i, j] <- sum(mat[row_indices, col_indices])
#     }
#   }
#   
#   return(result)
# }
# 
# # Example usage
# # Create an example 837x837 matrix
# example_matrix <- matrix(1:(837*837), nrow = 837, ncol = 837)
# 
# # Define block size
# block_size <- 6
# 
# # Get the resulting coarser matrix
# coarser_matrix <- sum_blocks(example_matrix, block_size)
# 
# # Display the resulting matrix
# print(coarser_matrix)
