# Read data from file
require('epi')

fluData <- read.csv("data/call.csv", header = TRUE)
data <- fluData[,2]
# nSum <- 2
# data <- sumData(data, nSum)

# Fitting epidemics
startOffset <- 1
endOffset <- 1
offsets <- list(startOffset=startOffset, endOffset=endOffset, minTruncation=4)

# Threshold
thresholds <- list(lim=0.9)

# Init Params = beta, gamma, S0
initParams <- c();

epiTypes <- c(0)

# Init Conds = S0, I0, R0
# I0 from first data point
initConds <- c();
plotConfig <- list(title="Synthedemic Decomposition of Carly Rae Jepsen BitTorrent Downloads", fileName="output/graphs/callSingle/", dataFile="output/data/call/callSingle.RData", envFile="output/data/call/callEnv.RData", pat=12, rat=60)

# Fit parameters
fitOverTimeMulti("LMS", 1:length(data), data, initConds, initParams, epiTypes, offsets, thresholds, plotConfig)