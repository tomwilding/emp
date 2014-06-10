load("output/data/call/callData90S2.RData")
require("epi")

fluData <- read.csv("data/call.csv", header = TRUE)
data <- fluData[,2]

# nSum <- 2
# data <- sumData(data, nSum)
times <- c(1:length(data))

# Fitting epidemics
startOffset <- 1
endOffset <- 1
offsets <- list(startOffset=startOffset, endOffset=endOffset, minTruncation=4)

# Threshold
thresholds <- list(lim=0.9)

# Init Params = beta, gamma, S0
initParams <- c()
# Epidemic type array epidemic types correspond to the number of parameters of the sub epidemic model
epiTypes <- c(0)
# Init Conds = S0, I0, R0
# I0 from first data point
initConds <- c()

plotConfig <- list(title="Synthedemic Decomposition of Caly Rae Jepsen BitTorrent Downloads", fileName="output/graphs/callFinal290/", dataFile="output/data/call/callData.RData", envFile="output/data/call/callEnv.RData", pat=80, rat=500)

# Fit parameters
# reconstructPlot(times, data, offsets, thresholds, initParams, initConds, plotConfig)

# Plot RSq graph
# plotRSq(times, data, offsets, thresholds, initParams, initConds, plotConfig)

# Plot Residuals
# plotResiduals(times, data, offsets, thresholds, initParams, initConds, plotConfig)

# t+1 prediction
# plotPred(times, data, offsets, thresholds, initParams, initConds, plotConfig)
analysis(times, data, offsets, thresholds, initParams, initConds, plotConfig)