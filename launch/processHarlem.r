# Read data from file
require('epi')

fluData <- read.csv("data/bauer.csv", header = TRUE)
# nSum <- 2
data <- fluData[,2]
# data <- sumData(data, nSum)

# Fitting epidemics
startOffset <- 4
endOffset <- 1
offsets <- list(startOffset=startOffset, endOffset=endOffset, minTruncation=5)

# Threshold
thresholds <- list(lim=0.95)

# Init Params = beta, gamma, S0
initParams <- c()
# Epidemic type array epidemic types correspond to the number of parameters of the sub epidemic model
epiTypes <- c(0)
# Init Conds = S0, I0, R0
# I0 from first data point
initConds <- c()

plotConfig <- list(title="Synthedemic Decomposition of Harlem Shake YouTube Views", fileName="output/graphs/harlem2/", dataFile="output/data/harlem/harlemData.RData", envFile="output/data/harlem/callEnv.RData", pat=12, rat=60)

# Fit parameters
fitOverTimeMulti("LMS", 1:length(data), data, initConds, initParams, epiTypes, offsets, thresholds, plotConfig)