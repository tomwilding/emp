require('epi')
# Simulate data 
fluData <- simSIR(0.002,0.1,500,10)
# fluData1 <- simExp(0.2,400)
fluData1 <- simSIR(0.002,0.1,400,10)
# Get data from dataframe
# Ensure first is larger than second
# nSum <- 4
# positiveInfectious <- fluData$data[,3]
# positiveInfectious <- sumData(fluData$data[,3], nSum)
positiveInfectious1 <- fluData1$data[,3]
# positiveInfectious1 <- sumData(fluData1$data[,3], nSum)

# Offset of t0 for second epidemic
# offset1 <- 50
offset1 <- 60
# Total length of the combined data 
totalLength <- length(positiveInfectious1) + offset1
positiveInfectious <- array(0,totalLength)
# Padding of zeros to offset data
positiveInfectiousPad <- numeric(offset1)
# Combine data with padding offset zeros
allPositiveInfectious <- c(positiveInfectious,positiveInfectiousPad)
allPositiveInfectious1 <- c(positiveInfectiousPad,positiveInfectious1)
# Add together the different predicted infectious values truncated to required size
data <- (allPositiveInfectious[1:totalLength]) + (allPositiveInfectious1[1:totalLength])
# data <- positiveInfectious1
# data <- positiveInfectious
tmax <- length(data)
times <- c(1:tmax)
# print(maxt)

# Fitting epidemics
startOffset <- 1
endOffset <- 1
minTruncation <- 160
offsets <- list(startOffset=startOffset, endOffset=endOffset, minTruncation=minTruncation)

# Target rSquare error
target <- 0.9

# Init Params = beta, gamma, S0
# initParams <- c(log(0.001),log(0.1),log(data[startOffset]*10),0, log(0.001),log(0.1),log(data[startOffset]*10),log(60/108))
initParams <- c(log(0.001),log(0.1),log(data[offset1 + 1]*10),0)


# Epidemic type array epidemic types correspond to the number of parameters of the sub epidemic model
# epiTypes <- c(4, 4)
epiTypes <- c(0, 4)
# Init Conds = S0, I0, R0
# I0 from first data point
# initConds <- c(1,data[startOffset],0,0 ,1,1,0,0)
initConds <- c(1,data[offset1 + 1],0,0)

plotConfig <- list(title="Synthedemic Decomposition of Simulated Data", fileName="output/graphs/mix/", dataFile="output/data/mix/mixData.RData", envFile="output/data/mix/mixEnv.RData", pat=5, rat=30)

# Fit parameters
fitOverTimeMulti("LMS", times, data, initConds, initParams, epiTypes, offsets, target, plotConfig, tmax)