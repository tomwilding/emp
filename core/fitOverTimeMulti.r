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
	
	# Initialise other parameters
	rSquare <- 0
	totalRSquare <- 0
	
	# Step size for iterative fitting
	step <- 1
	# Initial t0 value
	ts <- c(1)#, 14, 72, 133, 203, 259)
	# ts <- c(1, 34, 93)
	# ts <- c(1, 13, 52)

	# Set the number of epidemics
	k <- length(ts)

	# All evaluation vector
	evalList <- c()
	
	# Number of increasing residuals
	nRes <- 2

	# Track epidemic start time
	# startTime <- ts[1]
	# startTimeCount <- 0
	timeSinceOutbreak <- 0

	# testParams(times, data, initConds, params, epiTypes, ts, k, 0.1)
	################################################# Decompose Epidemics ################################################
	# Truncate the data to i data points from 20 within offset data
	for (i in seq(from=minTruncation, to=maxTruncation, by=step)) {
		# Fit k epidemics
		print("------------------------------------------------", quote=FALSE)
		print(paste(c("fitting "," of "), c(i, maxTruncation)), quote=FALSE); print(paste("k", k), quote=FALSE)
		print(ts)
		# print(paste("Beta", exp(initParams[2])))
		epidemicType <- epiTypes[k]
		# Determine if current epidemic start time is set
		# startTimeCount <- countStartTime(ts[k], startTime, startTimeCount)
		startTime <- ts[k]
		# startTimePrev <- ts[max(1, (k - 1))]
		# print(paste("Count", startTimeCount))
		# Determine epidemic type and fit over required range
		# startSearch <- max(1, startTime - 10)
		# endSearch <- max(1, min((startTime + 10), (i - minTruncation)))
		# Spike Epidemic or No epidemic
		eval <- fitInRangeParallel(setSolver(optimMethod, k, epiTypes), i, offsetTimes, offsetData, initConds, initParams, epiTypes, ts, k, plotConfig, 1)
		# Update parameters
		# maxt <- eval$optimTime
		# ts[k] <- maxt
		rSquare <- eval$optimRSquare
		optimParams <- eval$optimParams
		optimConds <- eval$initConds
		print("optimParams")
		print(optimParams)
		print(paste("rs", rSquare))
		# Get last residual and update residuals vector
		residuals <- eval$residuals
		
		lim <- thresholds$lim
		print(paste("timeSinceOutbreak", timeSinceOutbreak))
		
		# Check for redundant epidemics
		if ((rSquare > lim) && (k > 2)) {
			prevEpidemicType <- epiTypes[k - 1]
			curEpidemicType <- epiTypes[k]
			print(">>> Fit k-1 epidemics", quote=FALSE)
			initParamsLess <- reduceParams(initParams, curEpidemicType)
			initCondsLess <- reduceParams(initConds, curEpidemicType)
			evalLess <- fitInRangeParallel(setSolver(optimMethod, k - 1, epiTypes[1:(k - 1)]), i, offsetTimes, offsetData, initCondsLess, initParamsLess, epiTypes[1:(k - 1)], ts[1:(k - 1)], k - 1, plotConfig, 0)		
			lessRSquare <- evalLess$optimRSquare
			print(paste("lrs", lessRSquare))
			print(paste("lim",lim))
			# outbreakInRange <- detectOutbreakInRange(minTruncation, evalLess$residuals, nRes, startTimePrev, k)
			# print(paste("OutbreakDetectedReduce", outbreakInRange))
			if (lessRSquare > lim) {
				print("reduce epidemics")
				# Set k - 1 epidemics from now on
				k <- k - 1
				# Update parameters to continue fitting with k - 1 epidemics
				ts <- ts[1 : k]
				rSquare <- lessRSquare
				initParams <- initParamsLess
				initConds <- initCondsLess
				epiTypes <- epiTypes[1 : k]
				eval <- evalLess				
			}
		}

		# Try to improve the fit if rSquare has deteriorated
		outbreak <- detectOutbreak(eval$residuals, nRes, startTime, k)
		print(paste("Outbreak", outbreak))
		if ((rSquare < lim) && (outbreak > 0)) {
			# Try k+1 epidemics
			print(">>> Fit k+1", quote=FALSE)
			if (outbreak == 4 || outbreak == 0) {
				# SIR Detected
				initParamsMore <- c(initParams, c(0,0,0,0))
				initCondsMore <- c(initConds, c(1,1,0,0))
				epiTypesMore <- c(epiTypes, 4)
				# Fit SIR epidemic starting minTruncation before detected time
				tsMore <- c(ts, i)
				evalMore <- fitInRangeParallel(setSolver(optimMethod, k + 1, epiTypesMore), i, offsetTimes, offsetData, initCondsMore, initParamsMore, epiTypesMore, tsMore, k + 1, plotConfig, 1)
				# TODO: Update all times from optimisation, not just last time
				RSquareMore <- evalMore$optimRSquare
			} else if (outbreak == 1) {
				# EXP Detected
				initParamsMore <- c(initParams, 0)
				initCondsMore <- c(initConds, 1)
				epiTypesMore <- c(epiTypes, 1)
				tsMore <- c(ts, i)
				# Fit more epidemics with t0 set at i
				evalMore <- fitInRangeParallel(setSolver(optimMethod, k + 1, epiTypesMore), i, offsetTimes, offsetData, initCondsMore, initParamsMore, epiTypesMore, tsMore, k + 1, plotConfig, 1)
				RSquareMore <- evalMore$optimRSquare
			}

			if (RSquareMore > rSquare) {	
				# Current epidemic start time is not set
				# startTimeCount <- 0
				timeSinceOutbreak <- 0
				# Set k+1 epidemics from now on
				k <- k + 1
				eval <- evalMore
				# Update parameters
				ts <- tsMore
				rSquare <- eval$optimRSquare
				initParams <- initParamsMore
				initConds <- initCondsMore
				epiTypes <- epiTypesMore
			}
			# prevParams <- optimParams
			# prevConds <- optimConds
		}
		# } else if ((rSquare < lim) && (epidemicType == 1)) {
		# 	print("Reset Exp start time")
		# 	# New spike outbreak detected - replace old I0
		# 	initConds[length(initConds)] <- offsetData[i]
		# 	ts[k] <- i
		# }

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


detectOutbreak <- function(residuals, nRes, startTime, k) {
	outbreak <- 0
	# Ensure more than one residual before the last n residuals to calculate sdRes
	resLength <- length(residuals)
	# Ensure time since outbreak is greater than minTruncation
	if (resLength > startTime + nRes + minTruncation) {
		# Get standard deviation of residuals before the ones considered
		inRangeResiduals <- residuals[startTime : (resLength - nRes)]
		# minRes <- max(1, resLength - window)
		# inRangeResiduals <- abs(residuals[(resLength - window):(resLength - nRes)])
		# print(inRangeResiduals)
		meanRes <- mean(inRangeResiduals)
		sdRes <- sd(inRangeResiduals)

		print(paste("meanRes", meanRes))
		print(paste("sdRes ", sdRes))
		print(paste("Outbreaklim", meanRes + sdRes * 2))
		print(paste("Explim", meanRes + sdRes * 6))

		# Check if last n residuals are above set number of sd
		# Index of first residual to check
		startResIndex <- resLength - nRes + 1
		# Take last residuals
		outbreakRes <- min(residuals[startResIndex : resLength])
		expRes <- residuals[resLength]
		print(paste("OutbreakRes", outbreakRes))
		print(paste("ExpRes", expRes))
		# Set epidemic type according to residual limit
		outbreakLim <- (meanRes + (sdRes * 2))
		expLim <- (meanRes + (sdRes * 6))
		# If minimum residual increase is more than required, then set type
		if (outbreakRes > outbreakLim) {
			outbreak <- 4
		} else if (expRes > expLim) {
			outbreak <- 1
		}
	}
	outbreak	
}

detectOutbreakInRange <- function(minTruncation, residuals, nRes, startTime, k) {
	outbreakInRange <- 0
	for (i in (startTime + minTruncation):length(residuals)) {
		outbreakInRange <- outbreakInRange || detectOutbreak(residuals[1:i], nRes, startTime, k)
	}
	outbreakInRange
}

reduceParams <- function(initVec, epidemicType) {
	newVec <- initVec[1:(length(initVec) - epidemicType)]
}

# countStartTime <- function(startTimeNew, startTime, startTimeCount) {
# 	if (startTimeNew == startTime) {
# 		startTimeCount <- startTimeCount + 1
# 	}
# 	startTimeCount
# }