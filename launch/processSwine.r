# Read data from file
# require('epi')

epiData <- read.csv("data/h1n1.csv", header = TRUE)
# nSum <- 4
data <- epiData[,2]

times <- c(1:length(data))
# plot(times, data)
# readline()
# Fitting epidemics
minTruncation <- 3
# startOffset <- findStartOffset(data, minTruncation)
startOffset <- 1
endOffset <- 1
offsets <- list(startOffset=startOffset, endOffset=endOffset, minTruncation=minTruncation)
# Threshold
thresholds <- list(lim=0.9)

# Init Params = beta, gamma, S0
initParams <- c(log(1))
# Epidemic type array epidemic types correspond to the number of parameters of the sub epidemic model
epiTypes <- c(1)

# Init Conds = S0, I0, R0
# I0 from first data point
initConds <- c(1)

plotConfig <- list(title="2009 H1N1 Influenza Outbreak", fileName="output/graphs/swineBaseline/", dataFile="output/data/swineBaseline.RData", envFile="output/data/blur2/blurEnv.RData", pat=5, rat=25)

# gradientSearch(times, data, plotConfig)
# readline()

# Fit parameters
fitOverTimeMulti("LMS", times, data, initConds, initParams, epiTypes, offsets, thresholds, plotConfig)