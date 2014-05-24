# Read data from file
# require('epi')

epiData <- read.csv("data/blurred_lines_1.csv", header = TRUE)
# nSum <- 4
data <- epiData[,2]
# data <- sumData(data, nSum)
times <- c(1:length(data))

# Fitting epidemics
minTruncation <- 4
# startOffset <- findStartOffset(data, minTruncation)
startOffset <- 1
endOffset <- 1
offsets <- list(startOffset=startOffset, endOffset=endOffset, minTruncation=minTruncation)

# Threshold
thresholds <- list(lim=0.9)

# Init Params = beta, gamma, S0
initParams <- c()
# initParams <- c(logit(1e-5, 1e-3, 1e-2), logit(1e-3, 1e-1, 0.5), logit(1e2, 1e3, 1e6), 0, 
# 				logit(1e-5, 1e-3, 1e-2), logit(1e-3, 1e-1, 0.5), logit(1e2, 1e3, 1e6), 0,
# 				logit(1e-3, 1e-2, 0.5),
# 				logit(1e-5, 1e-3, 1e-2), logit(1e-3, 1e-1, 0.5), logit(1e2, 1e3, 1e6), 0,
				# logit(1e-6, 1e-3, 1e-2), logit(1e-3, 0.01, 1e-1), logit(1e2, 1e3, 1e6), 0)
# Epidemic type array epidemic types correspond to the number of parameters of the sub epidemic model
epiTypes <- c(0)
# epiTypes <- c(0, 4, 4, 1, 4, 4)

# Init Conds = S0, I0, R0
# I0 from first data point
initConds <- c()
# initConds <- c(	1,1,0,0, 
# 				1,1,0,0,
# 				1,
# 				1,1,0,0,
				# 1,1,0,0)

plotConfig <- list(title="Synthedemic Decomposition of Robin Thicke BitTorrent Downloads", fileName="output/graphs/blurOT/", dataFile="output/data/blurt2/logisInit1.RData", envFile="output/data/blur2/blurEnv.RData", pat=12, rat=60)

# gradientSearch(times, data, plotConfig)
# readline()

# Fit parameters
fitOverTimeMulti("LMS", times, data, initConds, initParams, epiTypes, offsets, thresholds, plotConfig)