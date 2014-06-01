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
	timeSinceOutbreak <- 0

	################################################# Decompose Epidemics ################################################
	foreach (i in seq(from=minTruncation, to=maxTruncation, by=step)) %dopar%{
		# Fit k epidemics
		print("------------------------------------------------", quote=FALSE)
		print(paste(c("fitting "," of "), c(i, maxTruncation)), quote=FALSE); print(paste("k", k), quote=FALSE)
		print(ts)
		epidemicType <- epiTypes[k]
		startTime <- ts[k]
		optimParams <- initParams
		# for (o in 1:5) {
			# print(paste("optim", o))
			eval <- fitInRangeParallel(setSolver(optimMethod, k, epiTypes), i, offsetTimes, offsetData, initConds, optimParams, epiTypes, ts, k, plotConfig, 1)
		# optimParams <- eval$optimParams
		# }
		# Update parameters
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
			optimParamsLess <- initParamsLess
			evalLess <- fitInRangeParallel(setSolver(optimMethod, k - 1, epiTypes[1:(k - 1)]), i, offsetTimes, offsetData, initCondsLess, optimParamsLess, epiTypes[1:(k - 1)], ts[1:(k - 1)], k - 1, plotConfig, 0)		
			optimParamsLess <- evalLess$optimParams
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
				initCondsMore <- c(initConds, c(1,1,0,0))
				epiTypesMore <- c(epiTypes, 4)
				tsMore <- c(ts, i)
				# initParamsMore <- getInitParams(setSolver(optimMethod, k + 1, epiTypesMore), i, offsetTimes, offsetData, initCondsMore, initParams, epiTypesMore, tsMore, k + 1, plotConfig, data)
				initParamsMore <- c(initParams, c(log(0.001), log(0.01), log(1000), logit((i - 40), (i - 10), i)))
			} else if (outbreak == 1) {
				# EXP Detected
				initCondsMore <- c(initConds, 1)
				epiTypesMore <- c(epiTypes, 1)
				tsMore <- c(ts, i)
				initParamsMore <- c(initParams, log(0.01))
			}

			# Optimise K + 1 epidemics
			optimParamsMore <- initParamsMore
			evalMore <- fitInRangeParallel(setSolver(optimMethod, k + 1, epiTypesMore), i, offsetTimes, offsetData, initCondsMore, optimParamsMore, epiTypesMore, tsMore, k + 1, plotConfig, 1)
			optimParamsMore <- evalMore$optimParams

			if (evalMore$optimRSquare > rSquare) {	
				# Current epidemic start time is not set
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
		}

		# Save all evaluations
		evalList[[i]] <- eval
		totalRSquare <- totalRSquare + rSquare
		# Update the initial parameters for the next fitting
		timeSinceOutbreak <- timeSinceOutbreak + 1
	}
	
	# Save all params
	save(evalList, file=plotConfig$dataFile)
	avRS <- totalRSquare / length(seq(from=minTruncation, to=maxTruncation, by=step))
	print(avRS)
	avRS
}


detectOutbreak <- function(residuals, nRes, startTime, k) {
	outbreak <- 0
	# Ensure more than one residual before the last n residuals to calculate sdRes
	resLength <- length(residuals)
	# Ensure time since outbreak is greater than minTruncation
	if (resLength > startTime + nRes + minTruncation) {
		# Get standard deviation of residuals before the ones considered
		inRangeResiduals <- residuals[startTime : (resLength - nRes)]
		meanRes <- mean(inRangeResiduals)
		sdRes <- sd(inRangeResiduals)

		print(paste("meanRes", meanRes))
		print(paste("sdRes ", sdRes))
		print(paste("Outbreaklim", meanRes + sdRes * 2))
		print(paste("Explim", meanRes + sdRes * 6))

		# Check if last n residuals are above set number of standard deviation
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