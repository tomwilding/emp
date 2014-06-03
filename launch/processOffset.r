require('epi')
# Simulate data 
fluData1 <- simSIR(0.002,0.1,400,1)
# fluData1 <- simExp(0.2,400)
# fluData2 <- simSIR(0.002,0.1,500,10)
# Get data from dataframe
# Ensure first is larger than second
# nSum <- 4
positiveInfectious1 <- fluData1$data[,3]
# positiveInfectious <- sumData(fluData$data[,3], nSum)
# positiveInfectious2 <- fluData2$data[,3]
# positiveInfectious1 <- sumData(fluData1$data[,3], nSum)

# Offset of t0 for second epidemic
# offset1 <- 50
offset1 <- 50
# offset2 <- 100
# Total length of the combined data assuming last epidemic is longest time
totalLength <- offset1 + length(positiveInfectious1)
# Padding of zeros to offset data
positiveInfectiousPadStart1 <- array(0, offset1)
# Combine data with padding offset zeros
# Add together the different predicted infectious values truncated to required size
data <- c(positiveInfectiousPadStart1, positiveInfectious1)
# data <- positiveInfectious1
tmax <- length(data)
times <- c(1:tmax)

# Fitting epidemics
startOffset <- 1
endOffset <- 1
minTruncation <- 4
offsets <- list(startOffset=startOffset, endOffset=endOffset, minTruncation=minTruncation)

# Target rSquare error
thresholds <- list(lim=0.9)

# Epidemic type array epidemic types correspond to the number of parameters of the sub epidemic model
epiTypes <- c(0)
# epiTypes <- c(4)

initParams <- c()
# initParams <- c(log(0.001), log(0.1), log(10), 0)
# Init Conds = S0, I0, R0
# I0 from first data point
initConds <- c()
# initConds <- c(1,log(data[startOffset + 1]),0,0)

plotConfig <- list(title="Synthedemic Decomposition of Simulated Data", fileName="output/graphs/offset/", dataFile="output/data/offset/offsetData.RData", envFile="output/data/offset/offsetEnv.RData", pat=5, rat=30)
# Fit parameters
fitOverTimeMulti("LMS", times, data, initConds, initParams, epiTypes, offsets, thresholds, plotConfig)