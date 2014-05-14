# require('epi')
# Simulate data 
fluData <- simSIR(0.001,0.05,400,1)
# fluData1 <- simExp(0.2,200)
fluData1 <- simSIR(0.001,0.1,400,1)
# Get data from dataframe
# Ensure first is larger than second
# nSum <- 4
positiveInfectious <- fluData$data[,3]
positiveInfectious1 <- fluData1$data[,3]
# positiveInfectious1 <- fluData1$data[,2]

# Offset of t0 for second epidemic
offset <- 30
offset1 <-50

# Padding of zeros to offset data
set.seed(1)
positiveInfectiousPadStart <- runif(offset)*2
positiveInfectiousPad1Start <- numeric(offset + offset1)
positiveInfectious1 <- c(positiveInfectiousPad1Start,positiveInfectious1)
totalLength <- length(positiveInfectious1)
# Combine data with padding offset zeros
positiveInfectious <- c(positiveInfectiousPadStart, positiveInfectious)
positiveInfectiousPadEnd <- c()
if (totalLength - length(positiveInfectious) > 0) {	
	positiveInfectiousPadEnd <- numeric(totalLength - length(positiveInfectious))
}
positiveInfectious <- c(positiveInfectious, positiveInfectiousPadEnd)

# Add together the different predicted infectious values truncated to required size
data <- positiveInfectious + positiveInfectious1

# data <- allPositiveInfectious

# Fitting epidemics
startOffset <- 1
endOffset <- 1
minTruncation <- 4
offsets <- list(startOffset=startOffset, endOffset=endOffset, minTruncation=minTruncation)

# Thresholds
thresholds <- list(lim=0.998)

# Init Params = beta, gamma, S0
initParams <- c();
# Epidemic type array epidemic types correspond to the number of parameters of the sub epidemic model
epiTypes <- c(0)
# Init Conds = S0, I0, R0
# I0 from first data point
initConds <- c();

plotConfig <- list(title="Synthedemic Decomposition of Simulated Data", fileName="output/graphs/mix/", dataFile="output/data/mix/mixData.RData", envFile="output/data/mix/mixEnv.RData", pat=5, rat=30)

# Fit parameters
fitOverTimeMulti("LMS", c(1:length(data)), data, initConds, initParams, epiTypes, offsets, thresholds, plotConfig)