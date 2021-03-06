load("output/data/harlem/harlemData.RData")
require("epi")

# Read data from file
epiData <- read.csv("data/harlem_shake.csv", header = TRUE)

# nSum <- 4
data <- epiData[,2]
# data <- sumData(data, nSum)
times <- c(1:length(data))

# Only fit over a specific range
startOffset <- 115
endOffset <- 1
offsets <- list(startOffset=startOffset, endOffset=endOffset, minTruncation=5, minTRange=3, maxTRange=3)

thresholds <- list(lim=0.9)

# Init Params = beta, gamma, S0
initParams <- c(log(0.001), log(0.1), log(data[startOffset]*10));

# Init Conds = S0, I0, R0
# I0 from first data point
initConds <- c(1,data[startOffset],0);

plotConfig <- list(title="Synthedemic Decomposition of Harlem Shake BitTorrent Downloads", fileName="output/graphs/harlem/", dataFile="output/data/harlem/harlemData.RData", envFile="output/data/harlem/blurEnv.RData", pat=12, rat=250)

# Fit parameters
# reconstructPlot(times, data, offsets, thresholds, initParams, initConds, plotConfig)

# Plot RSq graph
# plotRSq(times, data, offsets, thresholds, initParams, initConds, plotConfig)

# Plot Residuals
# plotResiduals(times, data, offsets, thresholds, initParams, initConds, plotConfig)

# t+1 prediction
plotPred(times, data, offsets, thresholds, initParams, initConds, plotConfig)