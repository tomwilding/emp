require(deSolve)

# Read in data
fluData <- read.csv("data/FluView_Data.csv", header = FALSE)

# Get flu data
percentPositive <- fluData[,11] / 100
totalSpecimins <- fluData[,3]
data <- totalSpecimins * percentPositive
# data <- c(numeric(30), data)
# Fitting epidemics
startOffset <- 7
endOffset <- 1
minTruncation <- 5
offsets <- list(startOffset=startOffset, endOffset=endOffset, minTruncation=minTruncation)
# Thresholds
thresholds <- list(diff=0.05, lim=0.9)
plotConfig <- list(title="Synthedemic Decomposition of CDC Data", fileName="output/graphs/CDCOT/", dataFile="output/data/CDCOT/CDCData.RData", envFile="output/data/CDC/CDCEnv.RData", pat=5, rat=30)

# Init Params = beta, gamma, S0
# Maybe update in iteration to start at the previous optim params in the next iteration?
# Assume at first point there are significantly more Susceptible than Infected?
initParams <- c(initParams <- c(logit(1e-5, 5e-4, 1e-2), logit(1e-3, 1e-1, 0.5), logit(1e2, incOrderOf(data[startOffset]), 1e6), 0))

epiTypes <- c(4)

# Init Conds = S0, I0, R0
# I0 from first data point
initConds <- c(1,data[startOffset],0,0);

# Fit parameters at different times
fitOverTimeMulti("LMS", c(1:length(data)), data, initConds, initParams, epiTypes, offsets, thresholds, plotConfig)