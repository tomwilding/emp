# Read data from file
require('epi')

fluData <- read.csv("data/harlem_shake_1.csv", header = TRUE)
# nSum <- 2
data <- fluData[,2]

# Fitting epidemics
startOffset <- 4
endOffset <- 1
offsets <- list(startOffset=startOffset, endOffset=endOffset, minTruncation=5)

# Threshold
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

plotConfig <- list(title="Synthedemic Decomposition of Harlem Shake YouTube Views", fileName="output/graphs/harlem1/", dataFile="output/data/harlem/harlemData.RData", envFile="output/data/harlem/callEnv.RData", pat=12, rat=60)

# Fit parameters
fitOverTimeMulti("LMS", 1:length(data), data, initConds, initParams, epiTypes, offsets, thresholds, plotConfig)