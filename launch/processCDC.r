require(deSolve)

# Read in data
fluData <- read.csv("data/FluView_Data.csv", header = FALSE)

# Get flu data
percentPositive <- fluData[,11] / 100
totalSpecimins <- fluData[,3]
data <- totalSpecimins * percentPositive

# Fitting epidemics
startOffset <- 7
endOffset <- 12
minTruncation <- 3
offsets <- list(startOffset=startOffset, endOffset=endOffset, minTruncation=minTruncation)
# Thresholds
thresholds <- list(diff=0.05, lim=0.9)
plotConfig <- list(title="Synthedemic Decomposition of CDC Data", fileName="output/graphs/CDC/", dataFile="output/data/CDC/CDCData.RData", envFile="output/data/CDC/CDCEnv.RData", pat=5, rat=30)

# Init Params = beta, gamma, S0
# Maybe update in iteration to start at the previous optim params in the next iteration?
# Assume at first point there are significantly more Susceptible than Infected?
initParams <- c(log(0.001), log(0.1), log(data[startOffset]*10));

# Init Conds = S0, I0, R0
# I0 from first data point
initConds <- c(1,data[startOffset],0);

# Fit parameters at different times
fitOverTimeMulti("LMS", c(1:length(data)), data, initConds, initParams, offsets, thresholds, plotConfig)