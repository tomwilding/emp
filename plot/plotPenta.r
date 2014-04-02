load("output/data/penta/pentaData.RData")
require("epi")

# Read data from file
fluData <- read.csv("data/penta.csv", header = TRUE)
# data <- sumData(fluData[,2], 4)
data <- fluData[,2]
times <- c(1:length(data))
# print(data)
# readline()
# times <- takeEveryOther(fluData[,1])
# times <- takeEveryOther(times)

# Fitting multiple epidemics
# load("data/pentaData.RData")
# Only fit over a specific range
startOffset <- 7
endOffset <- 1
offsets <- list(startOffset=startOffset, endOffset=endOffset, minTruncation=10, minTRange=3, maxTRange=3)

# Threshold
# thresholds <- list(diff=0.05, lim=0.96)
# thresholds <- list(diff=0.05, lim=0.7)
thresholds <- list(lim=0.9)
# Init Params = beta, gamma, S0
initParams <- c(log(0.005), log(0.5), log(data[startOffset]*10));

# Init Conds = S0, I0, R0
# I0 from first data point
initConds <- c(1,data[startOffset],0);

plotConfig <- list(title="Synthedemic Decomposition of Pentatonix Downloads", fileName="output/graphs/penta/", dataFile="output/data/penta/pentaData.RData", envFile="output/data/penta/pentaEnv.RData", pat=12, rat=250)

# Fit parameters
reconstructPlot(times, data, offsets, thresholds, initParams, initConds, plotConfig)

# Plot RSq graph
# plotRSq(times, data, offsets, thresholds, initParams, initConds, plotConfig)

# Plot Residuals
# plotResiduals(times, data, offsets, thresholds, initParams, initConds, plotConfig)

# t+1 prediction
# plotPred(times, data, offsets, thresholds, initParams, initConds, plotConfig)