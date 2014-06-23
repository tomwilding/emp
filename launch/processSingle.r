require('epi')
# Simulate data 
epiData <- simSIR(0.001,0.1,500,1)

# Get data from dataframe
positiveInfectious <- epiData$data[,3]

# Offset of t0 for second epidemic
offset <- 30

# Padding of zeros to offset data
positiveInfectiousPadStart <- numeric(offset)

# positiveInfectious <- c(positiveInfectiousPadStart, positiveInfectious)

# Add together the different predicted infectious values truncated to required size
data <- positiveInfectious

# data <- allPositiveInfectious

# Fitting epidemics
startOffset <- 1
endOffset <- 1
minTruncation <- 4
offsets <- list(startOffset=startOffset, endOffset=endOffset, minTruncation=minTruncation)

# Thresholds
thresholds <- list(lim=0.95)

# Init Params = beta, gamma, S0
# initParams <- c();
initParams <- c(log(0.001), log(0.01), log(1000))
# Epidemic type array epidemic types correspond to the number of parameters of the sub epidemic model
epiTypes <- c(3)
# Init Conds = S0, I0, R0
# I0 from first data point
# initConds <- c();
initConds <- c(1,1,0)

plotConfig <- list(title="Synthedemic Decomposition of Simulated Single Epidemic Data", fileName="output/graphs/single/", dataFile="output/data/single/singleData.RData", envFile="output/data/mix/mixEnv.RData", pat=5, rat=30)

# Fit parameters
fitOverTimeMulti("LMS", c(1:length(data)), data, initConds, initParams, epiTypes, offsets, thresholds, plotConfig)