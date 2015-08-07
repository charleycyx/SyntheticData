# -----------------------------------------------------
# Define Function: mtcdForEntireDataFrame
# -----------------------------------------------------

# This function simply calls the function above, ensuring that all the data is
# synthesized by setting the upperlimit equal to the highest count value

# -----------------------------------------------------
# 1. Roxygen Description
# -----------------------------------------------------

#' Create Synthetic Dataset
#'
#' Returns synthetic datasets from a Bayesian hierarchical model. This function wraps C++
#' code proposed in "Synthesizing Truncated Count Data for Confidentiality," and developed
#' by Sam Hawala, Jerry Reiter and Quanli Wang. References can be found in the package
#' description.
#' 
#'@param df A dataframe with three columns, the last column containing count to be synthesized
#'@param seed positive integer, random seed
#'@param niters positive integer, number of iterations for parameter generation
#'@param stride positive integer
#'@param upperLimit a positive defining the "small count" threshold
#'@param numModel a positive integer specifying the number of synthetic datasets to be generated
#'@param burnin a positive integer specifying the number of burnins
#'@param verbose logical. Set = TRUE for output debug information
#'
#'@section Details: [Add model detail here]
#'
#'@examples
#' getSyntheticData(df)
#' getSyntheticData(df, upperLimit = 10)
#' 
#' @return A dataframe containing the synthetic datasets and disclosure risk measures

# -----------------------------------------------------
# 2. Function Code
# -----------------------------------------------------

mtcdForEntireDataFrame <- function(df, seed=0, niters=100000, burnin=30000, stride=500, numModel=10, verbose=FALSE) {
  
  upperLimit = max(df[[3]])
  
  getSyntheticData(df,seed,niters,burnin,stride,numModel,verbose,upperLimit)
  
}
