# Read data from file
require('epi')

fluData <- read.csv("data/blurred_lines.csv", header = TRUE)
# nSum <- 4
data <- fluData[,2]
times <- (1:length(data))
# data <- sumData(data, nSum)

# Fitting epidemics
minTruncation <- 290
# startOffset <- findStartOffset(data, minTruncation)
startOffset <- 1
endOffset <- 1
offsets <- list(startOffset=startOffset, endOffset=endOffset, minTruncation=minTruncation)

# Threshold
thresholds <- list(lim=0.9)

# Init Params = beta, gamma, S0
initParams <- c(log(0.001),log(0.1),log(10),log(1),-3, 
				log(0.001),log(0.1),log(10),log(1),-0.5,
				log(0.001),log(0.1),log(10),log(1),1.5)
# initParams <- c(log(0.001),log(0.1),log(10),log(1),0)
# initParams <- c()

# Init Conds = S0, I0, R0
initConds <- c(1,1,0,0,0 ,1,1,0,0,0, 1,1,0,0,0)
# initConds <- c(1,1,0,0,0)
# initConds <- c()

# Epidemic type array epidemic types correspond to the number of parameters of the sub epidemic model
epiTypes <- c(0, 5, 5, 5)
# epiTypes <- c(0, 5)
# epiTypes <- c(0)

plotConfig <- list(title="Synthedemic Decomposition of Robin Thicke BitTorrent Downloads", fileName="output/graphs/blurOT/", dataFile="output/data/blurOT/blurData.RData", envFile="output/data/blur2/blurEnv.RData", pat=12, rat=60)
gradientSearch(times, data)
readline()
# Fit parameters
fitOverTimeMulti("LMS", times, data, initConds, initParams, epiTypes, offsets, thresholds, plotConfig)