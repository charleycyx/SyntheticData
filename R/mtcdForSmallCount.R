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
# 1. Function definition: mtcdForSmallCount
# -----------------------------------------------------

mtcdForSmallCount <- function(df, seed=0, niters=100000, burnin=30000, stride=500, numModel=10, verbose=FALSE, upperLimit=9) {
  
  # Create synthetic data for small counts
  dfs <- subset(df,df[[3]]<=upperLimit) 
  dfsyn <- getSyntheticData(dfs,seed,niters,burnin,stride,numModel,verbose,upperLimit)
  
  # Extract large count data
  dfcom <- subset(df,df[[3]]>upperLimit)
  
  #fill in all the synthetic data for not synthetic columns, using the original data
  for (i in 1:numModel) {
    dfcom <- cbind(dfcom,dfcom[[3]])
  }
  
  #fill in all the rest of the columns with NA
  for (i in 1:(length(dfsyn)-length(dfcom))) {
    dfcom <- cbind(dfcom,rep(NA,length(dfcom[[1]])))
  }
  
  #sync the names and rbind two data frame, this get us a data frame of the same size of the original one
  names(dfcom) <- names(dfsyn)
  dfsyn <- rbind(dfsyn,dfcom)
  
  #sync the names in df with the new df, so that we can merge for the original sequence
  names(df) <- names(dfsyn)[1:length(df)]
  
  merge(df[,1:2],dfsyn,sort=FALSE,all.x=TRUE,by=names(dfsyn)[1:2])
}
