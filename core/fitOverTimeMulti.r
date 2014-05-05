fitOverTimeMulti <- function(optimMethod, times, data, initConds, initParams, epiTypes, offsets, thresholds, plotConfig) {
	
	# Unpack starting parameters, conditions and offsets
	startOffset <- offsets$startOffset
	endOffset <- offsets$endOffset

	# Take data set within specified offset
	offsetTimes <- times[startOffset:(length(times)-endOffset)]
	offsetData <- data[startOffset:(length(data)-endOffset)]

	# Max and min truncated data set sizes within offset data
	minTruncation <- offsets$minTruncation
	maxTruncation <- length(offsetData)
	
	# Min and max t0 values to explore within offset data
	minTRange <- 3
	maxTRange <- 3
	
	# Initialise other parameters
	rSquare <- 0
	rSquareRefit <- 0 
	totalRSquare <- 0
	# Step size for iterative fitting
	step <- 1
	# Initial t0 value
	ts <- c(1)
	# Set the number of epidemics
	k <- 1
	# Number of data points to before epidemic
	window <- 20

	# All evaluation vector
	evalList <- c()
	
	# Number of increasing residuals
	nRes <- 2
	residuals <- c()

	# Track epidemic start time
	startTimePrev <- ts[1]
	startTimeCount <- 10
	repeatStartTimes <- 5
	timeSinceOutbreak <- 0

	epiList <- numeric(minTruncation)
	
	################################################# Decompose Epidemics ################################################
	# Truncate the data to i data points from minTruncation within offset data
	for (i in seq(from=minTruncation, to=maxTruncation, by=step)) {
		# Fit k epidemics
		print(paste("k",k))
		print("### Fit k", quote=FALSE); print(paste(c("fitting "," of "), c(i, maxTruncation)), quote=FALSE)
		# print(paste("Beta", exp(initParams[2])))
		epidemicType <- epiTypes[k]
		# Determine if current epidemic start time is set
		startTimeCount <- countStartTime(ts[k], startTimePrev, startTimeCount)
		startTimePrev <- ts[k]
		print(paste("Count", startTimeCount))
		# Determine epidemic type and fit over required range
		if (epidemicType == 3) {
			# SIR Epidemic
			if ((startTimeCount < repeatStartTimes) && (k > 1)) {
				# Explore t0 from previous epidemic start point
				eval <- fitInRangeParallel(setSolver(optimMethod, k, epiTypes), i, offsetTimes, offsetData, initConds, initParams, epiTypes, ts, k, c(minTruncation:(i - minTruncation)), plotConfig, 1)
			} else {
				# Epidemic starts at this time point
				eval <- fitInRangeParallel(setSolver(optimMethod, k, epiTypes), i, offsetTimes, offsetData, initConds, initParams, epiTypes, ts, k, c(ts[k]:ts[k]), plotConfig, 1)
			}
		} else {
			# Spike Epidemic or No epidemic
			eval <- fitInRangeParallel(setSolver(optimMethod, k, epiTypes), i, offsetTimes, offsetData, initConds, initParams, epiTypes, ts, k, c(ts[k]:ts[k]), plotConfig, 1)
		}
		# Update parameters
		maxt <- eval$optimTime
		ts[k] <- maxt
		rSquare <- eval$optimRSquare
		optimParams <- eval$multiParams
		optimConds <- eval$initConds
		print("optimParams")
		print(optimParams)
		# Get last residual and update residuals vector
		residuals <- c(residuals, eval$finalResidual)
		
		lim <- thresholds$lim
		epidemicType <- getEpidemicType(residuals, nRes, rSquare, window, k)
		print(paste("timeSinceOutbreak", timeSinceOutbreak))
		
		# Check for redundant epidemics
		if (rSquare > lim && k > 1) {
			prevEpidemicType <- epiTypes[k - 1]
			curEpidemicType <- epiTypes[k]
			print(">>> Fit k-1 epidemics", quote=FALSE)
			initParamsLess <- reduceParams(initParams, curEpidemicType)
			initCondsLess <- reduceParams(initConds, curEpidemicType)
			evalLess <- fitInRangeParallel(setSolver(optimMethod, k - 1, epiTypes[1:(k - 1)]), i, offsetTimes, offsetData, initCondsLess, initParamsLess, epiTypes[1:(k - 1)], ts[1:(k - 1)], k - 1, (ts[k - 1]:ts[k - 1]), plotConfig, 0)		
			lessRSquare <- evalLess$optimRSquare
			print(paste("lrs", lessRSquare))
			print(lim)
			print(lessRSquare > lim)
			if (lessRSquare > rSquare) {
				print("reduce epidemics")
				startTimeCount <- 0
				timeSinceOutbreak <- 0
				# Set k+1 epidemics from now on
				k <- k - 1
				# Update parameters to continue fitting with k+1 epidemics
				maxt <- evalLess$optimTime
				ts <- ts[1 : k]
				rSquare <- lessRSquare
				initParams <- initParamsLess
				initConds <- initCondsLess
				epiTypes <- epiTypes[1 : k]
				eval <- evalLess				
			}
		} 

		# Try to improve the fit if rSquare has deteriorated
		if ((timeSinceOutbreak > minTruncation) && (rSquare > 0 && rSquare < lim) || (epidemicType > 0)) {
			if (epidemicType == 0) { 
				epidemicType <- 3
			}
			# Set new epidemic type in epidemic type array
			epiTypesMulti <- c(epiTypes, epidemicType)
			# Try k+1 epidemics
			print(">>> Fit k+1", quote=FALSE)
			print(paste("!!!Epidemic", epidemicType), quote=FALSE)
			# Set initial parameters
			if (k == 1) {
				startParams <- c(log(0.001), log(0.01), log(data[i]*10))
				startConds <- c(1,data[i],0)
			}
			initParamsMulti <- newParams(optimParams, startParams, epidemicType)
			initCondsMulti <- newParams(optimConds, startConds, epidemicType)
			print("k+1 params")
			print(initParamsMulti)
			# Explore start range of SIR
			if (epidemicType == 3) {
				# Fit k + 1 epidemics with new SIR sub epidemic exploring t0
				evalMulti <- fitInRangeParallel(setSolver(optimMethod, k + 1, epiTypesMulti), i, offsetTimes, offsetData, initCondsMulti, initParamsMulti, epiTypesMulti, ts, k + 1, c(minTruncation:(i - minTruncation)), plotConfig, 1)
			} else {
				# Fit k + 1 epidemics with set t0 at i
				evalMulti <- fitInRangeParallel(setSolver(optimMethod, k + 1, epiTypesMulti), i, offsetTimes, offsetData, initCondsMulti, initParamsMulti, epiTypesMulti, ts, k + 1, c(i:i), plotConfig, 1)
			}
			# TODO: Does the residuals array need to be updated here? = Only consider residuals in initial fitInRangeParallel above - will be updated in next loop
			# residuals <- nIncResiduals[(i-(nResiduals-1)):i]
			multiRSquare <- evalMulti$optimRSquare
			# If K + 1 epidemic is better then set k + 1 from now on
			print(multiRSquare)
			print(rSquare)
			if (multiRSquare > rSquare) {
			print("epiSet")
			# if (multiRSquare > rSquare + (1 - lim)) {
				# Current epidemic start time is not set
				startTimeCount <- 0
				timeSinceOutbreak <- 0
				# Set k+1 epidemics from now on
				k <- k + 1
				# Update parameters to continue fitting with k+1 epidemics
				maxt <- evalMulti$optimTime
				ts <- c(ts, maxt)
				rSquare <- multiRSquare
				initParams <- initParamsMulti
				initConds <- initCondsMulti
				epiTypes <- epiTypesMulti
				eval <- evalMulti
			}
		} else if ((rSquare < lim) && (epidemicType == 1)) {
			print("Reset Exp start time")
			# New spike outbreak detected - replace old I0
			initConds[length(initConds)] <- offsetData[i]
			ts[k] <- i
		}

		# TODO: Store eval starts at i index causing NA
		evalList[[i]] <- eval
		totalRSquare <- totalRSquare + rSquare
		# Update the initial parameters for the next fitting
		# initParams <- multiParams
		# Increment time since outbreak
		timeSinceOutbreak <- timeSinceOutbreak + 1
	}
	
	# Save all params
	save(evalList, file=plotConfig$dataFile)
	avRS <- totalRSquare / length(seq(from=minTruncation, to=maxTruncation, by=step))
	print(avRS)
	avRS
	# save.image(plotConfig$envFile)
}

# Compare the last n residuals to sd of previous residuals before them
getEpidemicType <- function(residuals, nRes, rSquare, window, k) {
	incRes <- c()
	type <- 0
	# Standard deviation of residuals
	sdRes <- 0
	# Get the last n residuals
	residuals <- (residuals[!is.na(residuals)])

	lowerLimit <- 10

	# Ensure more than one residual before the last n residuals to calculate sdRes
	resLength <- length(residuals)
	# if (resLength > nRes + 1) {
	if (resLength > window + 1) {
		# Get standard deviation of residuals before the ones considered
		# absResiduals <- abs(residuals[1:(resLength - nRes)])
		minRes <- max(1, resLength - window)
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
		print(paste("SIRlim", meanRes + sdRes * 3))
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
			sirRes <- min(residuals[startResIndex:resLength])
			finalRes <- residuals[resLength]
			print(paste("FRes", finalRes))
			print(paste("SIRes", sirRes))
			# print(paste("SameSign",sameSign))
			# Set epidemic type according to residual limit
			# sirSD <- meanRes + (sdRes * 2)
			# spikeSD <- meanRes + (sdRes * 4)
			spikeLim <- meanRes + sdRes * 100
			sirLim <- meanRes + sdRes * 3
			# If minimum residual increase is more than required, then set type
			if ((finalRes > spikeLim) && (finalRes > lowerLimit)) {
				print("Spike Triggered")
				type <- 1
			} else if ((sirRes > sirLim) && sameSign) {
				print("SIR Triggered")
				type <- 3
			}
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