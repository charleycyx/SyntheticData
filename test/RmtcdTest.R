library(Rmtcd)

#do this to set the working directory to the current directory
script.dir <- dirname(sys.frame(1)$ofile)
setwd(script.dir)

#read the first data frame, where all cells are smaller than 9
read.table("data.txt") -> df

#synthesize the whole frame
mtcdForEntireDataFrame(df) -> resultEntireDf

#synthesize for only the smaller counts (<=5)
mtcdForSmallCount(df,upperLimit = 5) -> resultSmallCountDf

#read the second data frame, which has cells with counts as large as 20
read.table("data2.txt") -> df2

#same things as above, small counts defined as <=9
mtcdForEntireDataFrame(df2) -> resultEntireDf2
mtcdForSmallCount(df2,upperLimit = 9) -> resultSmallCountDf2

#read a data frame with county names(which should be converted into factors by read.table)
read.table("dataNames.txt") -> dfC

#same things as above
mtcdForEntireDataFrame(dfC) -> resultEntireDfN
mtcdForSmallCount(dfC,upperLimit = 5) -> resultSmallCountDfN