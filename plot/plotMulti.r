load("output/data/sim/simData.RData")


# ################################## Simulate data ########################################
fluData <- sim(0.001,0.1,500,10)
fluData1 <- sim(0.002,0.2,600,10)
# Get data from dataframe
positiveInfectious <- fluData$data[,3]
positiveInfectious1 <- fluData1$data[,3]
# Find max length as sim may return different lengths
maxLength <- max(length(positiveInfectious), length(positiveInfectious1))
# Offset of t0 for second epidemic
offset1 <- 40
# Total length of the combined data
totalLength <- maxLength + offset1
# Build times array
# times <- c(1:(totalLength));
# Padding of zeros to offset data
positiveInfectiousPad1 <- numeric(offset1)
# positiveInfectiousPad2 <- numeric(offset2)
# Combine data with padding offset zeros
allPositiveInfectious <- c(positiveInfectious,positiveInfectiousPad1)
allPositiveInfectious1 <- c(positiveInfectiousPad1,positiveInfectious1)
# Add together the different predicted infectious values truncated to required size
data <- (allPositiveInfectious[1:totalLength]) + (allPositiveInfectious1[1:totalLength])
# data <- sumData(data, 2)
times <- 1:length(data)
# fluData <- sim(0.002,0.1,500,10)
# fluData1 <- sim(0.002,0.2,400,10)
# fluData2 <- sim(0.005,0.1,200,10)
# # Get data from dataframe
# positiveInfectious <- (fluData$data[,3])
# positiveInfectious1 <- (fluData1$data[,3])
# positiveInfectious2 <- (fluData2$data[,3])
# # Find max length as sim may return different lengths
# maxLength <- max(length(positiveInfectious), length(positiveInfectious1))
# maxLength <- max(maxLength, length(positiveInfectious2))
# # Offset of t0 for second epidemic
# # offset1 <- 20
# offset1 <- 15
# offset2 <- 35
# # Total length of the combined data
# totalLength <- maxLength + offset2
# # Build times array
# times <- c(1:(totalLength));
# # Padding of zeros to offset data
# positiveInfectiousPad1 <- numeric(offset1)
# positiveInfectiousPad2 <- numeric(offset2)
# # Combine data with padding offset zeros
# allPositiveInfectious <- c(positiveInfectious,positiveInfectiousPad1)
# allPositiveInfectious1 <- c(positiveInfectiousPad1,positiveInfectious1)
# allPositiveInfectious2 <- c(positiveInfectiousPad2,positiveInfectious2)
# # Add together the different predicted infectious values truncated to required size
# allData <- (allPositiveInfectious[1:totalLength]) + (allPositiveInfectious1[1:totalLength]) + (allPositiveInfectious2[1:totalLength]);
# data <- takeEveryOther(allData)


################################## Fitting multiple epidemics ###########################
# Only fit over a specific range of times startOffset>=1
startOffset <- 1
endOffset <- 1
minTruncation <- 6
offsets <- list(startOffset=startOffset, endOffset=endOffset, minTruncation=minTruncation)

# Thresholds
thresholds <- list(diff=0.05, lim=0.95)

# Init Params = beta, gamma, S0
initParams <- c(log(0.001), log(0.1), log(data[startOffset]*10))

# Init Conds = S0, I0, R0
# I0 from first data point
initConds <- c(1,data[startOffset],0);

plotConfig <- list(title="Synthedemic Decomposition of Simulated Data", fileName="output/graphs/sim/", dataFile="output/data/sim/simData.RData", pat=5, rat=30)

reconstructPlot(times, data, offset, thresholds, initParams, initConds, plotConfig)
