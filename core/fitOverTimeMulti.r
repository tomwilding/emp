fitOverTimeMulti <- function(optimMethod, times, data, initConds, initParams, epiTypes, offsets, target, plotConfig, tmax) {
	
	# Unpack starting parameters, conditions and offsets
	startOffset <- offsets$startOffset
	endOffset <- offsets$endOffset

	# Take data set within specified offset
	offsetTimes <- times[startOffset:(length(times)-endOffset)]
	offsetData <- data[startOffset:(length(data)-endOffset)]

	# Max and min truncated data set sizes within offset data
	minTruncation <- offsets$minTruncation
	maxTruncation <- length(offsetData)
	
	# Initial t0 value
	# ts <- c(1, 23)
	# Set the number of epidemics
	k <- 2
	
	evalList <- c()

	################################################# Decompose Epidemics ################################################
	# Truncate the data to i data points from minTruncation within offset data
	for (i in seq(from=minTruncation, to=maxTruncation, by=1)) {
		print("### Fit k", quote=FALSE); print(paste(c("fitting "," of "), c(i, maxTruncation)), quote=FALSE)
		
		# Optimise k epidemics
		eval <- fitInRangeParallel(setSolver(optimMethod, k, epiTypes), i, offsetTimes, offsetData, initConds, initParams, epiTypes, k, tmax)
		eval <- fitInRangeParallel(setSolver(optimMethod, k, epiTypes), i, offsetTimes, offsetData, initConds, eval$multiParams, epiTypes, k, tmax)
		# maxt <- eval$optimTime
		# ts[k] <- maxt
		rSquare <- eval$optimRSquare
		multiParams <- eval$multiParams
		print(multiParams)
		print(rSquare)
		# TODO: Store eval starts at i index causing NA
		evalList[[i]] <- eval

		# # Check for redundant epidemics

		# # If rSquare has deteriorated then fit k + 1 epidemics
		# if (rSquare < target) {
		# 	# Update k
		# 	k <- k + 1
			
		# 	# Update parameters

		# 	# Optimise k + 1 epidemics

		# }
		# Set start parameters to current optimised parameters
		# initParams <- multiParams
	}
	
	# Save all params
	save(evalList, file=plotConfig$dataFile)
	# save.image(plotConfig$envFile)
}

# Compare the last n residuals to sd of previous residuals before them
getEpidemicType <- function(residuals, nRes, window, rSquare) {
	incRes <- c()
	type <- 0
	# Standard deviation of residuals
	sdRes <- 0
	# Lower limit of residual to determine outbreak
	lowerLimit <- 1000
	# Get the last n residuals
	residuals <- (residuals[!is.na(residuals)])
	
	# Ensure more than one residual before the last n residuals to calculate sdRes
	resLength <- length(residuals)
	# if (resLength > nRes + 1) {
	if (resLength > window + 1) {
		# Get standard deviation of residuals before the ones considered
		# absResiduals <- abs(residuals[1:(resLength - nRes)])
		absResiduals <- abs(residuals[(resLength - window):(resLength - nRes)])
		meanRes <- myMean(absResiduals)
		sdRes <- mySd(absResiduals, meanRes)
		# Reset between epidemics
		# meanDiffRes <- myMeanDiff(absResiduals)
		# sdDiffRes <- mySdDiff(absResiduals, meanDiffRes)
		# diffRes <- abs(residuals[resLength - nRes + 1] - residuals[resLength - nRes])
		print(paste("RSq", rSquare))
		print(paste("meanRes", meanRes))
		print(paste("sdRes ", sdRes))
		# print(paste("MeanDiffRes", meanDiffRes))
		# print(paste("SdDiffRes", sdDiffRes))
		# print(paste("DiffRes", diffRes))
		# If current residual sd is above zero check if last n residuals are above set number of sd
		if (sdRes > 0) {
			# Index of first residual to check
			startResIndex <- resLength - nRes + 1
			# Assume incRes is True and check condition for all n residuals
			sameSign <- max(residuals[startResIndex:resLength]) < 0 || min(residuals[startResIndex:resLength]) > 0
			# continuouslyIncreasing <- sort(residuals[startResIndex:resLength]) == residuals[startResIndex:resLength]
			sirIncRes <- min(abs(residuals[startResIndex:resLength]))
			finalRes <- abs(residuals[resLength])
			print(finalRes)
		}
		# print(paste("SameSign",sameSign))
		# Set epidemic type according to residual limit
		# sirSD <- meanRes + (sdRes * 2)
		# spikeSD <- meanRes + (sdRes * 4)
		spikeSD <- meanRes + sdRes * 6
		sirSD <- meanRes + sdRes * 2
		# If minimum residual increase is more than required, then set type
		if (finalRes > spikeSD && finalRes > lowerLimit) {
			print("Spike Set")
			type <- 1
		} else if (sirIncRes > sirSD && sirIncRes > lowerLimit && sameSign) {
			print("SIR Set")
			type <- 3
		}
	}
	type
}

newParams <- function(initVec, startVec, epidemicType) {
	if (epidemicType == 1) {
		# Update params for spike epidemic used I0 and gamma from initial
		newVec <- c(initVec, startVec[2])
	} else if (epidemicType == 3) {
		# Update params for SIR epidemic
		newVec <- c(initVec, startVec)
	}
}

reduceParams <- function(initVec, epidemicType) {
	if (epidemicType == 1) {
		# Update params for spike epidemic used I0 and gamma from initial
		newVec <- initVec[1:(length(initVec) - 1)]
	} else if (epidemicType == 3) {
		# Update params for SIR epidemic
		newVec <- initVec[1:(length(initVec) - 3)]
	}
}

countStartTime <- function(startTime, startTimePrev, startTimeCount) {
	if (startTime == startTimePrev) {
		startTimeCount <- startTimeCount + 1
	}
	startTimeCount
}