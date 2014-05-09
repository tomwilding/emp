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

	# All evaluation vector
	evalList <- c()
	
	# Number of increasing residuals
	nRes <- 2

	# Track epidemic start time
	startTime <- ts[1]
	startTimeCount <- 10
	repeatStartTimes <- 5
	timeSinceOutbreak <- 0

	epiList <- numeric(minTruncation)
	
	################################################# Decompose Epidemics ################################################
	# Truncate the data to i data points from minTruncation within offset data
	for (i in seq(from=minTruncation, to=maxTruncation, by=step)) {
		# Fit k epidemics
		print("------------------------------------------------", quote=FALSE)
		print(paste(c("fitting "," of "), c(i, maxTruncation)), quote=FALSE); print(paste("k", k), quote=FALSE)
		print(ts)
		# print(paste("Beta", exp(initParams[2])))
		epidemicType <- epiTypes[k]
		# Determine if current epidemic start time is set
		startTimeCount <- countStartTime(ts[k], startTime, startTimeCount)
		startTime <- ts[k]
		startTimePrev <- ts[max(1, (k - 1))]
		print(paste("Count", startTimeCount))
		# Determine epidemic type and fit over required range
		if (epidemicType == 3) {
			startSearch <- max(1, (startTime - minTruncation))
			print(paste("k range", c(startSearch:(startTime + minTruncation))))
			# SIR Epidemic
			eval <- fitInRangeParallel(setSolver(optimMethod, k, epiTypes), i, offsetTimes, offsetData, initConds, initParams, epiTypes, ts, k, c(startSearch:(startTime + minTruncation)), plotConfig, 1)
		} else {
			# Spike Epidemic or No epidemic
			eval <- fitInRangeParallel(setSolver(optimMethod, k, epiTypes), i, offsetTimes, offsetData, initConds, initParams, epiTypes, ts, k, c(startTime:startTime), plotConfig, 1)
		}
		# Update parameters
		maxt <- eval$optimTime
		ts[k] <- maxt
		rSquare <- eval$optimRSquare
		optimParams <- eval$multiParams
		optimConds <- eval$initConds
		print("optimParams")
		print(optimParams)
		print(paste("rs", rSquare))
		# Get last residual and update residuals vector
		residuals <- eval$residuals
		
		lim <- thresholds$lim
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
			print(paste("lim",lim))
			outbreakDetected <- detectOutbreak(evalLess$residuals, nRes, startTimePrev)
			print(paste("OutbreakDetectedReduce", outbreakDetected))
			if (lessRSquare > lim && !outbreakDetected) {
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
		outbreakDetected <- detectOutbreak(eval$residuals, nRes, startTime)
		print(paste("OutbreakDetected", outbreakDetected))
		if (((rSquare > 0 && rSquare < lim) || outbreakDetected) && timeSinceOutbreak > minTruncation) {
			# Try k+1 epidemics
			print(">>> Fit k+1", quote=FALSE)
			# Try SIR
			if (k == 1) {
				startParamsSIR <- c(log(0.001), log(0.1), log(orderOf(data[i])*10))
				startCondsSIR <- c(1,orderOf(data[i]),0)
			}
			initParamsSIR <- c(optimParams, startParamsSIR)
			initCondsSIR <- c(optimConds, startCondsSIR)
			epiTypesSIR <- c(epiTypes, 3)
			# Fit SIR epidemic searching t0 range
			print(paste("k+1 range", c(startTime:(i - minTruncation))))
			evalSIR <- fitInRangeParallel(setSolver(optimMethod, k + 1, epiTypesSIR), i, offsetTimes, offsetData, initCondsSIR, initParamsSIR, epiTypesSIR, ts, k + 1, c(startTime:(i - minTruncation)), plotConfig, 1)
			SIRRSquare <- evalSIR$optimRSquare
			print(paste("SIRRS", SIRRSquare))

			# Try EXP
			if (k == 1) {
				startParamsEXP <- c(log(1e-6))
				startCondsEXP <- c(orderOf(data[i]))
			}
			initParamsEXP <- c(optimParams, startParamsEXP)
			initCondsEXP <- c(optimConds, startCondsEXP)
			epiTypesEXP <- c(epiTypes, 1)		
			# Fit EXP epidemic with t0 set at i
			evalEXP <- fitInRangeParallel(setSolver(optimMethod, k + 1, epiTypesEXP), i, offsetTimes, offsetData, initCondsEXP, initParamsEXP, epiTypesEXP, ts, k + 1, c(startTime:i), plotConfig, 1)
			EXPRSquare <- evalEXP$optimRSquare
			print(paste("EXPRS", EXPRSquare))

			maxRSquare <- max(SIRRSquare, EXPRSquare)

			# If K + 1 epidemic is better then set k + 1 from now on
			print(paste("RS", rSquare))

			if (maxRSquare > rSquare) {	
				# Current epidemic start time is not set
				startTimeCount <- 0
				timeSinceOutbreak <- 0
				# Set k+1 epidemics from now on
				k <- k + 1
				if (SIRRSquare > EXPRSquare) {
					print("SIREpiSet")
					eval <- evalSIR
				} else {
					print("EXPEpiSet")
					eval <- evalEXP
				}
				# Update parameters
				ts <- c(ts, eval$optimTime)
				rSquare <- eval$optimRSquare
				initParams <- eval$multiParams
				initConds <- eval$initConds
				epiTypes <- eval$epiTypes
			}
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


detectOutbreak <- function(residuals, nRes, startTimePrev) {
	outbreak <- FALSE
	# Ensure more than one residual before the last n residuals to calculate sdRes
	resLength <- length(residuals)
	if (resLength > nRes + 1) {
		# Get standard deviation of residuals before the ones considered
		absResiduals <- abs(residuals[startTimePrev : (resLength - nRes)])
		# minRes <- max(1, resLength - window)
		# absResiduals <- abs(residuals[(resLength - window):(resLength - nRes)])
		meanRes <- myMean(absResiduals)
		sdRes <- mySd(absResiduals, meanRes)

		print(paste("meanRes", meanRes))
		print(paste("sdRes ", sdRes))
		print(paste("Outbreaklim", meanRes + sdRes * 6))
		print(paste("Explim", meanRes + sdRes * 10))
		# print(paste("MeanDiffRes", meanDiffRes))
		# print(paste("SdDiffRes", sdDiffRes))
		# print(paste("DiffRes", diffRes))
		# If current residual sd is above zero check if last n residuals are above set number of sd
		if (sdRes > 0) {
			# Index of first residual to check
			startResIndex <- resLength - nRes + 1
			# Assume incRes is True and check condition for all n residuals
			# sameSign <- max(residuals[startResIndex : resLength]) < 0 || min(residuals[startResIndex : resLength]) > 0
			# continuouslyIncreasing <- sort(residuals[startResIndex:resLength]) == residuals[startResIndex:resLength]
			outbreakRes <- min(residuals[startResIndex : resLength])
			expRes <- residuals[resLength]
			print(paste("OutbreakRes", outbreakRes))
			print(paste("ExpRes", expRes))
			# Set epidemic type according to residual limit
			outbreakLim <- (meanRes + (sdRes * 6))
			expLim <- (meanRes + (sdRes * 10))
			# If minimum residual increase is more than required, then set type
			outbreak <- ((outbreakRes > outbreakLim) || (expRes > expLim))
		}
	}
	outbreak	
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

countStartTime <- function(startTimeNew, startTime, startTimeCount) {
	if (startTimeNew == startTime) {
		startTimeCount <- startTimeCount + 1
	}
	startTimeCount
}