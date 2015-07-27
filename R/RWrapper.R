getSyntheticData <- function(df, seed=0, niters=100000, burnin=30000, stride=500, numModel=10, verbose=FALSE, upperLimit=9) {
  vdata <- numeric()
  for (i in 1:length(df[[1]])) {
    for (j in 1:3){
      vdata <- c(vdata,df[[j]][i])
    }
  }
  n <- length(df[[1]])
  odf <- data.frame(cbind(df,getSynData(vdata, n, seed, niters, burnin, stride, numModel, verbose, upperLimit-1)))
  
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