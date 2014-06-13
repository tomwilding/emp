require('epi')
# Simulate data 
epiData <- read.csv("data/flu.csv", header = TRUE)
# print(epiData)
# Get Infected number
data <- epiData$I
# plot(1:length(data), data)
# readline()

# Fitting epidemics
startOffset <- 1
endOffset <- 0
minTruncation <- 4
offsets <- list(startOffset=startOffset, endOffset=endOffset, minTruncation=minTruncation)

# Thresholds
thresholds <- list(lim=0.95)

# Init Params = beta, gamma, S0
# initParams <- c();
initParams <- c()
# Epidemic type array epidemic types correspond to the number of parameters of the sub epidemic model
epiTypes <- c(0)
# Init Conds = S0, I0, R0
# I0 from first data point
# initConds <- c();
initConds <- c()

plotConfig <- list(title="Initial Influenza Fitting", fileName="output/graphs/flu1/", dataFile="output/data/flu/fluData.RData", envFile="output/data/mix/mixEnv.RData", pat=5, rat=30)

# Fit parameters
fitOverTimeMulti("LMS", c(1:length(data)), data, initConds, initParams, epiTypes, offsets, thresholds, plotConfig)