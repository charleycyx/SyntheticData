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
