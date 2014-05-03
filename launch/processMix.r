require('epi')
if (file.exists("optimParams.txt")){file.remove("optimParams.txt")}
# Simulate data 
fluData1 <- simSIR(0.002,0.1,400,10)
# fluData1 <- simExp(0.2,400)
fluData2 <- simSIR(0.002,0.1,500,10)
# Get data from dataframe
# Ensure first is larger than second
# nSum <- 4
positiveInfectious1 <- fluData1$data[,3]
# positiveInfectious <- sumData(fluData$data[,3], nSum)
positiveInfectious2 <- fluData2$data[,3]
# positiveInfectious1 <- sumData(fluData1$data[,3], nSum)

# Offset of t0 for second epidemic
# offset1 <- 50
offset1 <- 50
offset2 <- 100
# Total length of the combined data assuming last epidemic is longest time
totalLength <- offset2 + length(positiveInfectious2)
positiveInfectious0 <- array(0, totalLength)
# Padding of zeros to offset data
positiveInfectiousPadStart1 <- array(0, offset1)
positiveInfectiousPadEnd1 <- array(0, totalLength - length(positiveInfectious1) - offset1)
positiveInfectiousPadStart2 <- array(0, offset2)
# Combine data with padding offset zeros
# Add together the different predicted infectious values truncated to required size
positiveInfectious1 <- c(positiveInfectiousPadStart1, positiveInfectious1, positiveInfectiousPadEnd1)
positiveInfectious2 <- c(positiveInfectiousPadStart2, positiveInfectious2)
data <- positiveInfectious0 + positiveInfectious1 + positiveInfectious2

# data <- positiveInfectious1
# data <- positiveInfectious
tmax <- length(data)
times <- c(1:tmax)
# print(maxt)

# Fitting epidemics
startOffset <- 1
endOffset <- 1
minTruncation <- 210
offsets <- list(startOffset=startOffset, endOffset=endOffset, minTruncation=minTruncation)

# Target rSquare error
target <- 0.95

# Init Params = beta, gamma, S0
initParams <- c(log(0.001),log(0.1),log(data[offset1 + 1]*10),log(data[offset1 + 1]),logit(50,tmax), log(0.001),log(0.1),log(data[offset1 + 1]*10),log(data[offset1 + 1]),logit(100,tmax))
# initParams <- c(log(0.001),log(0.1),log(100),log(1),0)
# initParams <- c()

# Init Conds = S0, I0, R0
initConds <- c(1,1,0,0,0 ,1,1,0,0,0)
# initConds <- c(1,1,0,0,0)
# initConds <- c()

# Epidemic type array epidemic types correspond to the number of parameters of the sub epidemic model
epiTypes <- c(0, 5, 5)
# epiTypes <- c(0, 5)
# epiTypes <- c(0)

plotConfig <- list(title="Synthedemic Decomposition of Simulated Data", fileName="output/graphs/mix/", dataFile="output/data/mix/mixData.RData", envFile="output/data/mix/mixEnv.RData", pat=5, rat=30)

# Fit parameters
fitOverTimeMulti("LMS", times, data, initConds, initParams, epiTypes, offsets, target, plotConfig)