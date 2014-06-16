load("output/data/swineSingle.RData")
require("epi")

# Read data from file
epiData <- read.csv("data/h1n1.csv", header = TRUE)

# nSum <- 4
data <- epiData[,2]
# data <- sumData(data, nSum)
times <- c(1:length(data))

# Only fit over a specific range
startOffset <- 1
endOffset <- 1
offsets <- list(startOffset=startOffset, endOffset=endOffset, minTruncation=4, minTRange=3, maxTRange=3)

thresholds <- list(lim=0.8)

# Init Params = beta, gamma, S0
initParams <- c();

# Init Conds = S0, I0, R0
# I0 from first data point
initConds <- c();

plotConfig <- list(title="Synthedemic Decomposition of the 2009 H1N1 Outbreak in England", fileName="output/graphs/report/swine/", dataFile="output/data/blur2/blurData90I.RData", envFile="output/data/blur2/blurEnv.RData", pat=5, rat=25)

# Fit parameters
# reconstructPlot(times, data, offsets, thresholds, initParams, initConds, plotConfig)

# Plot RSq graph
# plotRSq(times, data, offsets, thresholds, initParams, initConds, plotConfig)

# Plot Residuals
# plotResiduals(times, data, offsets, thresholds, initParams, initConds, plotConfig)

# t+1 prediction
plotPred(times, data, offsets, thresholds, initParams, initConds, plotConfig)
# analysis(times, data, offsets, thresholds, initParams, initConds, plotConfig)