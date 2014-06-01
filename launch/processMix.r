# Simulate data
fluData <- simSIR(0.001,0.05,400,1)
fluData1 <- simSIR(0.002,0.1,300,1)

positiveInfectious <- fluData$data[,3]
positiveInfectious1 <- fluData1$data[,3]

# Offset of t0 for second epidemic
offset <- 30
offset1 <- 45

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
times <- c(1:length(data))

# Fitting epidemics
startOffset <- 1
endOffset <- 1
minTruncation <- 4
offsets <- list(startOffset=startOffset, endOffset=endOffset, minTruncation=minTruncation)

# Thresholds
thresholds <- list(lim=0.995)

# Init Params = beta, gamma, S0
initParams <- c()
# initParams <- c(log(0.001), log(0.01), log(1000), logit((34 - 20), (34 - 10), 34),
				# log(0.001), log(0.01), log(1000), logit((85 - 20), (85 - 10), 85))
# Epidemic type array epidemic types correspond to the number of parameters of the sub epidemic model
epiTypes <- c(0)
# epiTypes <- c(0, 4, 4)

# Init Conds = S0, I0, R0
# I0 from first data point
initConds <- c()
# initConds <- c(1,1,0,0, 1,1,0,0)

plotConfig <- list(title="Synthedemic Decomposition of Simulated Data", fileName="output/graphs/mixOT/", dataFile="output/data/mix/mixDataFTt10.RData", envFile="output/data/mix/mixEnv.RData", pat=5, rat=30)

# Fit parameters
fitOverTimeMulti("LMS", times, data, initConds, initParams, epiTypes, offsets, thresholds, plotConfig)