load("output/data/blur/blurData.RData")

# Read data from file
fluData <- read.csv("data/blurred_lines.csv", header = TRUE)
data <- sumData(fluData[,2], 4)
times <- c(1:length(data))
# print(data)
# readline()
# times <- takeEveryOther(fluData[,1])
# times <- takeEveryOther(times)

# Fitting multiple epidemics
# load("data/blurData.RData")
# Only fit over a specific range
startOffset <- 1
endOffset <- 1
offsets <- list(startOffset=startOffset, endOffset=endOffset, minTruncation=14, minTRange=3, maxTRange=3)

# Threshold
# thresholds <- list(diff=0.05, lim=0.96)
# thresholds <- list(diff=0.05, lim=0.7)
thresholds <- list(diff=0.12, lim=0.9)
# Init Params = beta, gamma, S0
initParams <- c(log(0.001), log(0.1), log(data[startOffset]*10));

# Init Conds = S0, I0, R0
# I0 from first data point
initConds <- c(1,data[startOffset],0);

plotConfig <- list(title="Synthedemic Decomposition of Robin Thicke BitTorrent Downloads", fileName="output/graphs/blur/", dataFile="output/data/blur/blurData.RData", envFile="output/data/blur/blurEnv.RData", pat=12, rat=60)

# Fit parameters
reconstructPlot(times, data, offset, thresholds, initParams, initConds, plotConfig)

# Plot RSq graph
plotRSq(times, data, offset, thresholds, initParams, initConds, plotConfig)