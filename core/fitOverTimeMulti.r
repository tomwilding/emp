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
	ts <- c(1)

	# Set the number of epidemics
	k <- length(ts)

	# All evaluation vector
	evalList <- c()
	
	# Number of increasing residuals
	nRes <- 2

	# Track epidemic start time
	startTime <- ts[1]
	startTimeCount <- 0
	timeSinceOutbreak <- 0

	# testParams(times, data, initConds, params, epiTypes, ts, k, granularity)
	
	################################################# Iterative Epidemic Fitting ################################################
	# Truncate the data to i data points from 20 within offset data
	for (i in seq(from=minTruncation, to=maxTruncation, by=step)) {
		# Fit k epidemics
		print("------------------------------------------------", quote=FALSE)
		print(paste(c("fitting "," of "), c(i, maxTruncation)), quote=FALSE); print("Number of Epidemics:", quote=FALSE); print("", quote="FALSE"); print(k)
		print("", quote="FALSE")
		print("Epidemic start times:", quote=FALSE)
		print(ts)
		print("", quote="FALSE")
		epidemicType <- epiTypes[k]
		# Determine if current epidemic start time is set
		startTimeCount <- countStartTime(ts[k], startTime, startTimeCount)
		startTime <- ts[k]
		startTimePrev <- ts[max(1, (k - 1))]
		# Determine epidemic type and fit over required range
		startSearch <- max(1, startTime - 20)
		endSearch <- max(1, min((startTime + 20), (i - minTruncation)))
		if (epidemicType == 3) {
			print(paste("k range:", startSearch, "to ", endSearch), quote=FALSE)
			# SIR Epidemic
			eval <- fitInRangeParallel(setSolver(optimMethod, k, epiTypes), i, offsetTimes, offsetData, initConds, initParams, epiTypes, ts, k, c(startSearch:endSearch), plotConfig, 1)
		} else {
			# Spike Epidemic or No epidemic
			eval <- fitInRangeParallel(setSolver(optimMethod, k, epiTypes), i, offsetTimes, offsetData, initConds, initParams, epiTypes, ts, k, c(startTime:startTime), plotConfig, 1)
		}
		# Update parameters
		maxt <- eval$optimTime
		ts[k] <- maxt
		rSquare <- eval$optimRSquare
		optimParams <- eval$multiParams
		optimConds <- eval$multiConds
		print("", quote="FALSE")
		print("optimParams:", quote=FALSE)
		print(exp(optimParams), quote=FALSE)
		print("", quote="FALSE")
		print(paste("rs:", rSquare), quote=FALSE)
		# Get last residual and update residuals vector
		residuals <- eval$residuals
		
		lim <- thresholds$lim
		print(paste("timeSinceOutbreak:", timeSinceOutbreak), quote=FALSE)
		
		# Check for redundant epidemics
		if ((rSquare > lim) && (k > 2)) {
			prevEpidemicType <- epiTypes[k - 1]
			curEpidemicType <- epiTypes[k]
			print(">>> Fit k-1 epidemics", quote=FALSE)
			initParamsLess <- reduceParams(initParams, curEpidemicType)
			initCondsLess <- reduceParams(initConds, curEpidemicType)
			evalLess <- fitInRangeParallel(setSolver(optimMethod, k - 1, epiTypes[1:(k - 1)]), i, offsetTimes, offsetData, initCondsLess, initParamsLess, epiTypes[1:(k - 1)], ts[1:(k - 1)], k - 1, c(ts[k - 1]:ts[k - 1]), plotConfig, 0)		
			lessRSquare <- evalLess$optimRSquare
			print(paste("lrs", lessRSquare))
			print(paste("lim",lim))

			if (lessRSquare > lim) {
				print("reduce epidemics")
				# Set k - 1 epidemics from now on
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
		outbreak <- detectOutbreak(eval$residuals, nRes, startTime, k)
		print(paste("Outbreak:", outbreak), quote=FALSE)
		if ((rSquare < lim) && (outbreak > 0) && (timeSinceOutbreak > minTruncation)) {
			# Try k+1 epidemics
			print(">>> Fit k+1", quote=FALSE)
			if (outbreak == 3 || outbreak == 0) {
				# Try SIR
				initParamsMore <- c(initParams, c(log(0.001), log(0.01), log(1000)))
				initCondsMore <- c(initConds, c(1,1,0))
				epiTypesMore <- c(epiTypes, 3)
				# Fit More epidemic searching t0 range from previous epidemic
				startSearch <- startTime
				endSearch <- max(1, (i - minTruncation))
				print(paste("k+1 range:", startSearch, "to ", endSearch), quote=FALSE)
				evalMore <- fitInRangeParallel(setSolver(optimMethod, k + 1, epiTypesMore), i, offsetTimes, offsetData, initCondsMore, initParamsMore, epiTypesMore, ts, k + 1, c(startSearch:endSearch), plotConfig, 1)
				RSquareMore <- evalMore$optimRSquare

			} else if (outbreak == 1) {
				# EXP Detected
				initParamsMore <- c(initParams, log(0.01))
				initCondsMore <- c(initConds, 1)
				epiTypesMore <- c(epiTypes, 1)		
				# Fit More epidemic with t0 set at i
				evalMore <- fitInRangeParallel(setSolver(optimMethod, k + 1, epiTypesMore), i, offsetTimes, offsetData, initCondsMore, initParamsMore, epiTypesMore, ts, k + 1, c(startTime:i), plotConfig, 1)
				RSquareMore <- evalMore$optimRSquare
			}

			if (RSquareMore > rSquare) {	
				# Current epidemic start time is not set
				startTimeCount <- 0
				# Update time since outbreak for SIR
				timeSinceOutbreak <- 0
				# Set k+1 epidemics from now on
				k <- k + 1
				eval <- evalMore
				# Update parameters
				ts <- c(ts, eval$optimTime)
				rSquare <- eval$optimRSquare
				initParams <- initParamsMore
				initConds <- initCondsMore
				epiTypes <- epiTypesMore
			}
		}

		evalList[[i]] <- eval
		totalRSquare <- totalRSquare + rSquare
		# Increment time since outbreak
		timeSinceOutbreak <- timeSinceOutbreak + 1
	}
	
	# Save all params
	save(evalList, file=plotConfig$dataFile)
	avRS <- totalRSquare / length(seq(from=minTruncation, to=maxTruncation, by=step))
	print(avRS)
	avRS
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

countStartTime <- function(startTimeNew, startTime, startTimeCount) {
	if (startTimeNew == startTime) {
		startTimeCount <- startTimeCount + 1
	}
	startTimeCount
}