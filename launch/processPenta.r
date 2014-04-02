# Read data from file
require('epi')

fluData <- read.csv("data/penta.csv", header = TRUE)
# nSum <- 4
data <- fluData[,2]

# Fitting epidemics

startOffset <- 7
endOffset <- 1
offsets <- list(startOffset=startOffset, endOffset=endOffset, minTruncation=10)

# Threshold
# thresholds <- list(diff=0.05, lim=0.96)
# thresholds <- list(diff=0.05, lim=0.7)
thresholds <- list(lim=0.9)

# Init Params = beta, gamma, S0
initParams <- c(log(0.005), log(0.5), log(data[startOffset]*10))
# Epidemic type array epidemic types correspond to the number of parameters of the sub epidemic model
epiTypes <- c(3)
# Init Conds = S0, I0, R0
# I0 from first data point
initConds <- c(1,data[startOffset],0)

plotConfig <- list(title="Synthedemic Decomposition of Pentatonix Downloads", fileName="output/graphs/penta/", dataFile="output/data/penta/pentaData.RData", envFile="output/data/penta/pentaEnv.RData", pat=12, rat=60)

# Fit parameters
fitOverTimeMulti("LMS", 1:length(data), data, initConds, initParams, epiTypes, offsets, thresholds, plotConfig)