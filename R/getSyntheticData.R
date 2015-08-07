# -----------------------------------------------------
# Define Function: getSyntheticData
# -----------------------------------------------------

# This function does all the work and serves as a base function which is called
# by the other exported functions. After reformatting the input dataframe,
# it calls the C++ model (getSynData) and then formats the output.

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
# 1. Function definition: getSyntheticData
# -----------------------------------------------------


getSyntheticData <- function(df, seed=0, niters=100000, burnin=30000, stride=500, numModel=10, verbose=FALSE, upperLimit=9) {
  
  # Create vector (vdata) by layering each row of the 3-col input dataframe alongside each other
  n <- length(df[[1]])
  vdata <- numeric()				
  for (i in 1:n) {
    for (j in 1:3){
      vdata <- c(vdata,df[[j]][i])
    }
  }

  # Run the dataframe through the C++ model and append the results to the input dataframe
  odf <- data.frame(cbind(df,getSynData(vdata, n, seed, niters, burnin, stride, numModel, verbose, upperLimit-1)))
  
  # Name the the dataframe variables
  synNames <- character(numModel)
  for (i in 1:numModel) {
    synNames[i] <- paste("syntheticData",as.character(i),sep="")
  }
  riskNames <- character(numModel*3)
  for (i in 0:(numModel-1)) {
    riskNames[3*i+1] <- paste("dt",as.character(i),sep="")
    riskNames[3*i+2] <- paste("Rall",as.character(i),sep="")
    riskNames[3*i+3] <- paste("Runq",as.character(i),sep="")
  }
  
  names(odf) <- c("source", "destination", "y", synNames, "CIlower","CIupper",riskNames)
  odf
}
