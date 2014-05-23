# ################################## Simulate data ########################################
fluData <- simSIR(0.001,0.1,500,1)
fluData1 <- simSIR(0.002,0.2,600,1)
# Get data from dataframe
positiveInfectious <- fluData$data[,3]
positiveInfectious1 <- fluData1$data[,3]

# Offset of t0 for second epidemic
offset <- 10
offset1 <- 35

# Build times array
# Padding of zeros to offset data
positiveInfectiousPad1 <- numeric(offset1)

# Combine data with padding offset zeros
allPositiveInfectious1 <- c(positiveInfectiousPad1,positiveInfectious1)

offsetEnd <- length(allPositiveInfectious1) - length(positiveInfectious)
positiveInfectiousPadEnd <- numeric(offsetEnd)
allPositiveInfectious <- c(positiveInfectious,positiveInfectiousPadEnd)
# Add together the different predicted infectious values truncated to required size
set.seed(1)
data <- (allPositiveInfectious) + (allPositiveInfectious1)
data <- c(runif(offset)*2,  data)
times <- c(1:length(data))

# Fitting epidemics
startOffset <- 1
endOffset <- 0
minTruncation <- 4
offsets <- list(startOffset=startOffset, endOffset=endOffset, minTruncation=minTruncation)

# Thresholds
thresholds <- list(lim=0.995)

# Init Params = beta, gamma, S0
initParams <- c()
# initParams <- c(0, 0, 0, 0,
# 				0, 0, 0, 0)
# Epidemic type array epidemic types correspond to the number of parameters of the sub epidemic model
epiTypes <- c(0)
# epiTypes <- c(0, 4, 4)

# Init Conds = S0, I0, R0
# I0 from first data point
initConds <- c()
# initConds <- c(1,1,0,0, 1,1,0,0);

plotConfig <- list(title="Synthedemic Decomposition of Simulated Data", fileName="output/graphs/mixOT1/", dataFile="output/data/mix/mixData1.RData", envFile="output/data/mix/mixEnv.RData", pat=5, rat=30)

# gradientSearch(times, data, plotConfig)
# readline()
# Fit parameters
fitOverTimeMulti("LMS", times, data, initConds, initParams, epiTypes, offsets, thresholds, plotConfig)