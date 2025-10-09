########################################################################################################
# It is a common issue in count matrix analysis specially for miRNAs to see duplicated gene names.
# So, you can set the gene names as the row names of the matrix and you can do some downstream analysis. 
# There are two options to deal with redundant genes in a count matrix: 
#      1) keep one gene with average count values. 
#      2) keep all genes by changing their names (for example by adding "_i" to the (i-1)th redundant gene)
# I used this simple function to find and correct the redundant gene names based on option 2.
#########################################################################################################

suffix_duplicates <- function(x) {
  # Check if the input is a character vector
  if (!is.character(x)) {
    stop("Input must be a character vector.")
  }

  # 'ave' is used to perform a cumulative count for each unique element in the vector.
  # For each group of identical elements (defined by 'x' itself), it generates a sequence.
  # For A = c("a","b","b","c","c","c"), this produces c(1, 1, 2, 1, 2, 3)
  counts <- ave(rep(1, length(x)), x, FUN = cumsum)

  # Create a logical vector to identify elements that are duplicates (i.e., count > 1)
  is_duplicate <- counts > 1

  # For the elements that are duplicates, create the suffix "_[count - 1]"
  # We subtract 1 because the first duplicate should have suffix "_1" (e.g., the 2nd 'b').
  # paste0 is a more efficient version of paste for concatenating strings.
  x[is_duplicate] <- paste0(x[is_duplicate], "_", counts[is_duplicate] - 1)

  # Return the modified vector
  return(x)
}
