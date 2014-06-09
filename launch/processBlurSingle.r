# Read data from file
require('epi')

fluData <- read.csv("data/blurred_lines.csv", header = TRUE)
nSum <- 4
data <- fluData[,2]
# data <- sumData(data, nSum)

# Fitting epidemics
startOffset <- 1
endOffset <- 1
offsets <- list(startOffset=startOffset, endOffset=endOffset, minTruncation=5)

# Threshold
thresholds <- list(lim=0.9)

# Init Params = beta, gamma, S0
initParams <- c()
# Epidemic type array epidemic types correspond to the number of parameters of the sub epidemic model
epiTypes <- c(0)
# Init Conds = S0, I0, R0
# I0 from first data point
initConds <- c()

plotConfig <- list(title="Synthedemic Decomposition of Robin Thicke BitTorrent Downloads", fileName="output/graphs/blurSingle/", dataFile="output/data/blurSingle/blurSingleData.RData", envFile="output/data/blurSingle/blurEnv.RData", pat=12, rat=60)

# Fit parameters
fitOverTimeSingle("LMS", 1:length(data), data, initConds, initParams, epiTypes, offsets, thresholds, plotConfig)