require('epi')
# Simulate data 
totalRSquare <- 0
i <- 0
for (bf in seq(from=0.001, to=0.006, by=0.001)) {
	for (gf in seq(from=0.1, to=0.3, by=0.05)) {
		for (bs in seq(from=0.001, to=0.006, by=0.001)) {
			for (gs in seq(from=0.1, to=0.3, by=0.05)) {
				fluData <- simSIR(bf, gf, 500, 10)
				fluData1 <- simSIR(bs, gs, 600, 10)

				# Get data from dataframe
				# Ensure first is larger than second
				nSum <- 4
				positiveInfectious <- sumData(fluData$data[,3], nSum)
				positiveInfectious1 <- sumData(fluData1$data[,3], nSum)
				# Offset of t0 for second epidemic
				offset1 <- 30
				# Total length of the combined data
				# Ensure length(positiveInfectious) < length(positiveInfectious1) + offset 
				totalLength <- length(positiveInfectious1) + offset1
				# Build times array
				times <- c(1:(totalLength));
				# Padding of zeros to offset data
				positiveInfectiousPad <- numeric(offset1)
				# Combine data with padding offset zeros
				allPositiveInfectious <- c(positiveInfectious,positiveInfectiousPad)
				allPositiveInfectious1 <- c(positiveInfectiousPad,positiveInfectious1)
				# Add together the different predicted infectious values truncated to required size
				data <- (allPositiveInfectious[1:totalLength]) + (allPositiveInfectious1[1:totalLength])
				print(bf);print(gf);print(bs);print(gs)

				# Fitting epidemics
				startOffset <- 1
				endOffset <- 1
				minTruncation <- 3
				offsets <- list(startOffset=startOffset, endOffset=endOffset, minTruncation=minTruncation)

				# Thresholds
				thresholds <- list(diff=0.05, lim=0.9)

				# Init Params = beta, gamma, S0
				initParams <- c(log(0.001), log(0.1), log(data[startOffset]*10));

				# Init Conds = S0, I0, R0
				# I0 from first data point
				initConds <- c(1,data[startOffset],0);

				plotConfig <- list(title="Synthedemic Decomposition of Simulated Data", fileName="output/graphs/sim/", dataFile=paste("output/data/sim/simData", bf, bs, bs, gs, ".RData"), envFile="output/data/sim/simEnv.RData", pat=5, rat=30, bf=bf, gf=gf, bs=bs, gs=gs)

				# Fit parameters
				totalRSquare <- totalRSquare + fitOverTimeMulti("LMS", c(1:length(data)), data, initConds, initParams, offsets, thresholds, plotConfig)
				i <- i + 1
				print(paste("################## Run ", i, " of 1296"))
			}
		}
	}
}
avRSquare <- totalRSquare / i
print(paste("Average RSquare: ", avRSquare, " over ", i, " runs"))

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