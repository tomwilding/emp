load("output/data/mix/mixDataFT140.RData")
require('epi')
# Simulate data 
# require('epi')
# Simulate data 
fluData <- simSIR(0.0001,0.05,4000,1)
fluData1 <- simSIR(0.0001,0.1,4000,1)
# Get data from dataframe
positiveInfectious <- fluData$data[,3]
positiveInfectious1 <- fluData1$data[,3]

# Offset of t0 for second epidemic
offset <- 30
offset1 <- 50

# Build times array
# Padding of zeros to offset data
# positiveInfectiousPad <- numeric(offset)
positiveInfectiousPad1 <- numeric(offset1)

# Combine data with padding offset zeros
allPositiveInfectious1 <- c(positiveInfectiousPad1, positiveInfectious1)


positiveInfectiousPadEnd <- numeric(length(allPositiveInfectious1))
allPositiveInfectious <- c(positiveInfectious, positiveInfectiousPadEnd)
# Add together the different predicted infectious values truncated to required size
totalLength <- length(allPositiveInfectious1)
set.seed(1)
data <- allPositiveInfectious + allPositiveInfectious1
# Offset start
data <- c(runif(offset)*2,  data[1:totalLength])
times <- c(1:length(data))

# Fitting epidemics
startOffset <- 1
endOffset <- 0
minTruncation <- 4
offsets <- list(startOffset=startOffset, endOffset=endOffset, minTruncation=minTruncation)

# Thresholds
thresholds <- list(lim=0.995)

# Init Params = beta, gamma, S0
initParams <- c()
# initParams <- c(log(0.001), log(0.01), log(1000),
# 				log(0.001), log(0.01), log(1000))
# Epidemic type array epidemic types correspond to the number of parameters of the sub epidemic model
epiTypes <- c(0)
# epiTypes <- c(0, 3, 3)

# Init Conds = S0, I0, R0
# I0 from first data point
initConds <- c()
# initConds <- c(1,1,0, 1,1,0);
plotConfig <- list(title="Synthedemic Decomposition of Simulated Data", fileName="output/graphs/mixOT5/", dataFile="output/data/mix/mixData.RData", envFile="output/data/mix/mixEnv.RData", pat=25, rat=150)

# Fit parameters
reconstructPlot(times, data, offsets, thresholds, initParams, initConds, plotConfig)

# Plot Residuals
# plotResiduals(times, data, offsets, thresholds, initParams, initConds, plotConfig)
# plotRSq(times, data, offsets, thresholds, initParams, initConds, plotConfig)

# t+1 prediction
# plotPred(times, data, offsets, thresholds, initParams, initConds, plotConfig)
# analysis(times, data, offsets, thresholds, initParams, initConds, plotConfig)