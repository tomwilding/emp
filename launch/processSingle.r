require('epi')
# Simulate data 
fluData <- simSIR(0.001,0.1,500,10)

# Get data from dataframe
positiveInfectious <- fluData$data[,3]

# Offset of t0 for second epidemic
offset <- 30

# Padding of zeros to offset data
positiveInfectiousPadStart <- numeric(offset)

positiveInfectious <- c(positiveInfectiousPadStart, positiveInfectious)

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
initParams <- c();
# Epidemic type array epidemic types correspond to the number of parameters of the sub epidemic model
epiTypes <- c(0)
# Init Conds = S0, I0, R0
# I0 from first data point
initConds <- c();

plotConfig <- list(title="Synthedemic Decomposition of Simulated Data", fileName="output/graphs/single/", dataFile="output/data/single/singleData.RData", envFile="output/data/mix/mixEnv.RData", pat=5, rat=30)

# Fit parameters
fitOverTimeMulti("LMS", c(1:length(data)), data, initConds, initParams, epiTypes, offsets, thresholds, plotConfig)