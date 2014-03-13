# Read data from file
require('epi')

fluData <- read.csv("data/blurred_lines.csv", header = TRUE)
# nSum <- 4
data <- fluData[,2]

# Fitting epidemics

startOffset <- 1
endOffset <- 1
offsets <- list(startOffset=startOffset, endOffset=endOffset, minTruncation=5)

# Threshold
# thresholds <- list(diff=0.05, lim=0.96)
# thresholds <- list(diff=0.05, lim=0.7)
thresholds <- list(lim=0.85)

# Init Params = beta, gamma, S0
initParams <- c(log(0.005), log(0.5), log(data[startOffset]*10));
# Epidemic type array epidemic types correspond to the number of parameters of the sub epidemic model
epiTypes <- c(3)
# Init Conds = S0, I0, R0
# I0 from first data point
initConds <- c(1,data[startOffset],0);

plotConfig <- list(title="Synthedemic Decomposition of Robin Thicke BitTorrent Downloads", fileName="output/graphs/blur1/", dataFile="output/data/blur1/blurData.RData", envFile="output/data/blur1/blurEnv.RData", pat=12, rat=60)

# Fit parameters
fitOverTimeMulti("LMS", 1:length(data), data, initConds, initParams, epiTypes, offsets, thresholds, plotConfig)