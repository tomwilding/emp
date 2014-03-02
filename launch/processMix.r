require('epi')
# Simulate data 
fluData <- simSIR(0.001,0.05,500,10)
fluData1 <- simExp(0.2,400)
# Get data from dataframe
# Ensure first is larger than second
positiveInfectious <- fluData$data[,3]
positiveInfectious1 <- fluData1$data[,2]


# Offset of t0 for second epidemic
offset1 <- 23
# Total length of the combined data 
totalLength <- length(positiveInfectious1) + offset1
# Padding of zeros to offset data
positiveInfectiousPad <- numeric(offset1)
# Combine data with padding offset zeros
allPositiveInfectious <- c(positiveInfectious,positiveInfectiousPad)
allPositiveInfectious1 <- c(positiveInfectiousPad,positiveInfectious1)
# Add together the different predicted infectious values truncated to required size
data <- (allPositiveInfectious[1:totalLength]) + (allPositiveInfectious1[1:totalLength])

# Fitting epidemics
startOffset <- 1
endOffset <- 1
minTruncation <- 6
offsets <- list(startOffset=startOffset, endOffset=endOffset, minTruncation=minTruncation)

# Thresholds
thresholds <- list(diff=0.05, lim=0.9)

# Init Params = beta, gamma, S0
initParams <- c(log(0.001), log(0.1), log(data[startOffset]*10));

# Init Conds = S0, I0, R0
# I0 from first data point
initConds <- c(1,data[startOffset],0);

plotConfig <- list(title="Synthedemic Decomposition of Simulated Data", fileName="output/graphs/mix/", dataFile="output/data/mix/mixData.RData", envFile="output/data/mix/mixEnv.RData", pat=5, rat=30)

# Fit parameters
fitOverTimeMulti("LMS", c(1:length(data)), data, initConds, initParams, offsets, thresholds, plotConfig)