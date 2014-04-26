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
	totalRSqaure <- 0
	# Step size for iterative fitting
	step <- 1
	# Initial t0 value
	ts <- c(1)
	# Set the number of epidemics
	k <- 0

	# All evaluation vector
	evalList <- c()
	
	# Number of increasing residuals
	nRes <- 2
	window <- 20
	residuals <- c()

	# At the start the start time is set to 0
	startTimePrev <- ts[1]
	startTimeCount <- 10
	repeatStartTimes <- 5
	timeSinceOutbreak <- 0

	epiList <- numeric(minTruncation)
	
	################################################# Decompose Epidemics ################################################
	# Truncate the data to i data points from minTruncation within offset data
	for (i in seq(from=minTruncation, to=maxTruncation, by=step)) {
			print("### Fit k", quote=FALSE); print(paste(c("fitting "," of "), c(i, maxTruncation)), quote=FALSE)
		if (k < 1) {
			# Take mean as prediction
			eval <- fitInRangeParallel(setSolver(optimMethod, k, epiTypes), i, offsetTimes, offsetData, initConds, initParams, epiTypes, ts, k, ts, plotConfig, 1)
			# Update fit properties
			maxt <- eval$optimTime
			ts[1] <- maxt
		} else {
			# Fit k epidemics (k >= 1)
			# print(paste("Beta", exp(initParams[2])))
			epidemicType <- epiTypes[k]
			# Determine if current epidemic start time is set
			startTimeCount <- countStartTime(ts[k], startTimePrev, startTimeCount)
			startTimePrev <- ts[k]
			print(paste("Count", startTimeCount))
			# Determine epidemic type and fit over required range
			if (epidemicType == 3) {
				if (startTimeCount > repeatStartTimes || k == 1) {
					# If only 1 epidemic assume it starts at given time
					eval <- fitInRangeParallel(setSolver(optimMethod, k, epiTypes), i, offsetTimes, offsetData, initConds, initParams, epiTypes, ts, k, c(ts[k]:ts[k]), plotConfig, 1)
				} else {
					# print(paste("times", ts))
					print(paste("min", i - window))
					print(paste("max", i - maxTRange))
					# Explore t0 from previous epidemic start point
					eval <- fitInRangeParallel(setSolver(optimMethod, k, epiTypes), i, offsetTimes, offsetData, initConds, initParams, epiTypes, ts, k, c((i - window):(i - maxTRange)), plotConfig, 1)
				}
			} else if (epidemicType == 1) {
				eval <- fitInRangeParallel(setSolver(optimMethod, k, epiTypes), i, offsetTimes, offsetData, initConds, initParams, epiTypes, ts, k, c(ts[k]:ts[k]), plotConfig, 1)
			}
			# Update fit properties
			maxt <- eval$optimTime
			ts[k] <- maxt
			print(paste("ts", ts))
		}
		rSquare <- eval$optimRSquare
		multiParams <- eval$multiParams

		# Get last residual and update residuals vector
		residuals <- c(residuals, eval$finalResidual)
		# residuals <- eval$residuals
		# print(residuals)

		# Print next prediction
		# print(eval$nextPred)

		lim <- thresholds$lim
		epidemicType <- getEpidemicType(residuals, nRes, window, rSquare)

		# Check for redundant epidemics
		if (rSquare > lim && k > 1) {
			prevEpidemicType <- epiTypes[k - 1]
			curEpidemicType <- epiTypes[k]
			print(">>> Fit k-1 epidemics", quote=FALSE)
			initParamsLess <- reduceParams(initParams, curEpidemicType)
			initCondsLess <- reduceParams(initConds, curEpidemicType)
			evalLess <- fitInRangeParallel(setSolver(optimMethod, k - 1, epiTypes[1:(k - 1)]), i, offsetTimes, offsetData, initCondsLess, initParamsLess, epiTypes[1:(k - 1)], ts[1:(k - 1)], k - 1, ts[k - 1]:ts[k - 1], plotConfig, 0)		
			lessRSquare <- evalLess$optimRSquare
			print(paste("lrs", lessRSquare))
			print(lim)
			print(lessRSquare > lim)
			if (lessRSquare > lim) {
				print("reduce epidemics")
				startTimeCount <- 0
				timeSinceOutbreak <- 0
				# Set k+1 epidemics from now on
				k <- k - 1
				# Update parameters to continue fitting with k+1 epidemics
				multiParams <- evalLess$multiParams
				maxt <- evalLess$optimTime
				ts <- ts[1:(k-1)]
				rSquare <- lessRSquare
				initConds <- initCondsLess
				epiTypes <- epiTypes[1:(k - 1)]
				eval <- evalLess				
			}
		# } else if ((rSquare < lim) && (epidemicType > 0)) {
		} else if ((rSquare < lim && k > 1) || epidemicType > 0) {
		# if (epidemicType > 0) {
			# Set new epidemic type in epidemic type array
			epiTypesMulti <- c(epiTypes, epidemicType)
			# Try k+1 epidemics
			print(">>> Fit k+1", quote=FALSE)
			print(paste("!!!Epidemic", epidemicType), quote=FALSE)

			# Get rid of initial data with no epidemic
			if (k == 0) {
				offsetTimes <- offsetTimes[(maxt - minTRange):length(offsetTimes)]
				offsetData <- offsetData[(maxt - minTRange):length(offsetData)]
				# Set parameters
				startParams <- c(log(0.001), log(0.01), log(offsetData[1]*10))
				initParamsMulti <- startParams
				startConds <- c(1,offsetData[1],0)
				initCondsMulti <- startConds
			} else {
				initParamsMulti <- newParams(initParams, startParams, epidemicType)
				initCondsMulti <- newParams(initConds, startConds, epidemicType)
			}

			print(k)
			print(initParams)
			print(initConds)
			print(ts)
			print(epiTypesMulti)

			# Search start time of k + 1 epidemic for SIR and not initial
			if ((epidemicType == 3 || epidemicType == 0) && k > 0) {
				# Fit k+1 epidemics with new SIR sub epidemic exploring t0 from previous epidemic start point
				evalMulti <- fitInRangeParallel(setSolver(optimMethod, k + 1, epiTypesMulti), i, offsetTimes, offsetData, initCondsMulti, initParamsMulti, epiTypesMulti, ts, k + 1, c((i - window):(i - maxTRange)), plotConfig, 1)
			# If Exp or initial SIR then don't search start times
			} else if (epidemicType == 1 || k == 0) {
				# Fit k+1 epidemics with set t0 at i
				evalMulti <- fitInRangeParallel(setSolver(optimMethod, k + 1, epiTypesMulti), i, offsetTimes, offsetData, initCondsMulti, initParamsMulti, epiTypesMulti, ts, k + 1, c(i:i), plotConfig, 1)
			}
			# TODO: Does the residuals array need to be updated here? = Only consider residuals in initial fitInRangeParallel above - will be updated in next loop
			# residuals <- nIncResiduals[(i-(nResiduals-1)):i]
			multiRSquare <- evalMulti$optimRSquare
			# If K + 1 epidemic is better then set k + 1 from now on
			print(multiRSquare)
			print(rSquare)
			if (multiRSquare > rSquare) {
				print("Set k+1")
				# Current epidemic start time is not set
				startTimeCount <- 0
				timeSinceOutbreak <- 0
				# Set k+1 epidemics from now on
				k <- k + 1

				# Update parameters to continue fitting with k+1 epidemics
				multiParams <- evalMulti$multiParams
				maxt <- evalMulti$optimTime
				ts[k] <- maxt
				rSquare <- multiRSquare
				initConds <- initCondsMulti
				epiTypes <- epiTypesMulti
				eval <- evalMulti
			}
		} else if ((rSquare < lim) && (epidemicType == 1)) {
			# New spike outbreak detected - replace old I0
			initConds[length(initConds)] <- offsetData[i]
			ts[k] <- i
		}
		# TODO: Store eval starts at i offset
		evalList[[i]] <- eval
		totalRSqaure <- totalRSqaure + rSquare
		# Update the initial parameters for the next fitting
		initParams <- multiParams
		# Increment time since outbreak
		timeSinceOutbreak <- timeSinceOutbreak + 1
	}
	
	# Save all params
	save(evalList, file=plotConfig$dataFile)
	avRS <- totalRSqaure / length(seq(from=minTruncation, to=maxTruncation, by=step))
	print(avRS)
	avRS
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
	if (resLength > nRes + 1) {
	# if (resLength > window + 1) {
		# Get standard deviation of residuals before the ones considered
		obsResiduals <- residuals[1:(resLength - nRes)]
		# absResiduals <- abs(residuals[(resLength - window):(resLength - nRes)])
		meanRes <- myMean(obsResiduals)
		sdRes <- mySd(obsResiduals, meanRes)
		# Reset between epidemics
		# meanDiffRes <- myMeanDiff(obsResiduals)
		# sdDiffRes <- mySdDiff(obsResiduals, meanDiffRes)
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
			sirIncRes <- min(residuals[startResIndex:resLength])
			finalRes <- residuals[resLength]
			print(paste("finalRes", finalRes))
			print(paste("sirIncRes", sirIncRes))
		}
		# print(paste("SameSign",sameSign))
		# Set epidemic type according to residual limit
		# sirLim <- meanRes + (sdRes * 2)
		# spikeLim <- meanRes + (sdRes * 4)
		spikeLim <- meanRes + sdRes * 6
		sirLim <- meanRes + sdRes * 2
		# If minimum residual increase is more than required, then set type
		if (finalRes > spikeLim && finalRes > lowerLimit) {
			print("Spike Set")
			type <- 1
		} else if (sirIncRes > sirLim && sirIncRes > lowerLimit && sameSign) {
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
	} else if (epidemicType == 3 || epidemicType == 0) {
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