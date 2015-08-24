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
#'@param df A dataframe with three columns, the last column containing count to be synthesized, the first two must be numeric, numbers that uniquely identify counties
#'@param seed positive integer for random number generation
#'@param niters positive integer, number of iterations for parameter generation
#'@param stride the model will save a model(set of parameters) for every stride number of models
#'@param numModel a positive integer specifying the number of synthetic datasets to be generated
#'@param burnin a positive integer specifying the number of burnins(iterations of parameters will not be saved)
#'
#'@section Details: The first three columns of output are the original dataset. The next numModel columns are the synthetic datasets. The next two colunms give the 95 percent confidence intevals estimated using all saved models. The following columns give the dt, Rall and Runq risk measurements.
#'
#'@examples
#' setwd("00_pkg_src/Rmtcd/test")
#' read.table("data.txt") -> dataFrame
#' mtcdForEntireDataFrame(dataFrame, numModel = 5)
#' 
#' @return A dataframe containing the synthetic datasets and disclosure risk measures

# -----------------------------------------------------
# 2. Function Code
# -----------------------------------------------------

mtcdForEntireDataFrame <- function(df, seed=0, niters=100000, burnin=30000, stride=500, numModel=10) {
  
  upperLimit = max(df[[3]])
  
  getSyntheticData(df,seed,niters,burnin,stride,numModel,FALSE,upperLimit)
  
}
