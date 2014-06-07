# Read data from file
require('epi')

fluData <- read.csv("data/call.csv", header = TRUE)
# nSum <- 2
data <- fluData[,2]
# data <- sumData(data, nSum)

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
plotConfig <- list(title="Synthedemic Decomposition of Carly Rae Jepsen BitTorrent Downloads", fileName="output/graphs/call2/", dataFile="output/data/call/callDataNoOBLim.RData", envFile="output/data/call/callEnv.RData", pat=12, rat=60)

# Fit parameters
fitOverTimeMulti("LMS", 1:length(data), data, initConds, initParams, epiTypes, offsets, thresholds, plotConfig)