# Read data from file
require('epi')

fluData <- read.csv("data/blurred_lines.csv", header = TRUE)
# nSum <- 4
data <- fluData[,2]
# data <- sumData(data, nSum)

# Fitting epidemics
minTruncation <- 5
# startOffset <- findStartOffset(data, minTruncation)
startOffset <- 1
endOffset <- 1
offsets <- list(startOffset=startOffset, endOffset=endOffset, minTruncation=minTruncation)

# Threshold
thresholds <- list(lim=0.9)

# Init Params = beta, gamma, S0
initParams <- c(log(0.001), log(0.01), log(data[startOffset]*10))
# Epidemic type array epidemic types correspond to the number of parameters of the sub epidemic model
epiTypes <- c(3)
# Init Conds = S0, I0, R0
# I0 from first data point
initConds <- c(1,data[startOffset],0)

plotConfig <- list(title="Synthedemic Decomposition of Robin Thicke BitTorrent Downloads", fileName="output/graphs/blur2/", dataFile="output/data/blur2/blurData.RData", envFile="output/data/blur2/blurEnv.RData", pat=12, rat=60)

# Fit parameters
fitOverTimeMulti("LMS", 1:length(data), data, initConds, initParams, epiTypes, offsets, thresholds, plotConfig)