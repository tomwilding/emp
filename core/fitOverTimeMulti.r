fitOverTimeMulti <- function(optimMethod, times, data, initConds, initParams, epiTypes, offsets, thresholds, plotConfig) {
	
	# Unpack starting parameters, conditions and offsets
	startOffset <- offsets$startOffset
	endOffset <- offsets$endOffset

	startParams <- initParams
	startConds <- initConds
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
	k <- 1

	# All evaluation vector
	evalList <- c()
	
	# Number of increasing residuals
	nRes <- 2
	window <- 10
	residuals <- c()

	epiList <- numeric(minTruncation)
	
	################################################# Decompose Epidemics ################################################
	# Truncate the data to i data points from minTruncation within offset data
	for (i in seq(from=minTruncation, to=maxTruncation, by=step)) {
		# Fit k epidemics
		print("### Fit k", quote=FALSE); print(paste(c("fitting "," of "), c(i, maxTruncation)), quote=FALSE)
		# print(paste("Beta", exp(initParams[2])))
		epidemicType <- epiTypes[k]
		if (epidemicType == 3) {
			if (k == 1) {
				# If only 1 epidemic assume it starts at given time
				eval <- fitInRangeParallel(setSolver(optimMethod, k, epiTypes), i, offsetTimes, offsetData, initConds, initParams, epiTypes, ts, k, c(ts[k]:ts[k]), plotConfig)
			} else {
				eval <- fitInRangeParallel(setSolver(optimMethod, k, epiTypes), i, offsetTimes, offsetData, initConds, initParams, epiTypes, ts, k, c(minTRange:(i-maxTRange)), plotConfig)
			}
		} else if (epidemicType == 1) {
			eval <- fitInRangeParallel(setSolver(optimMethod, k, epiTypes), i, offsetTimes, offsetData, initConds, initParams, epiTypes, ts, k, c(ts[k]:ts[k]), plotConfig)
		}
		# TODO: Store eval starts at i offset
		evalList[[i]] <- eval

		maxt <- eval$optimTime
		ts[k] <- maxt
		rSquare <- eval$optimRSquare
		multiParams <- eval$multiParams
		# Get last residual and update residuals vector
		residuals <- c(residuals, eval$finalResidual)


		# Try to improve the fit if rSquare has deteriorated
		lim <- thresholds$lim
		epidemicType <- getEpidemicType(residuals, nRes, window, rSquare)
		# epidemicType <- getEpidemicType(residuals, nRes, window, rSquare)
		if ((rSquare < lim) && (epidemicType > 0)) {
		# if (epidemicType > 4) {
			# Set new epidemic type in epidemic type array
			epiTypesMulti <- c(epiTypes, epidemicType)
			# Try k+1 epidemics
			print(">>> Fit k+1", quote=FALSE)
			print(paste("!!!Epidemic", epidemicType), quote=FALSE)
			# Set initial parameters
			initParamsMulti <- newParams(initParams, startParams, epidemicType)
			initCondsMulti <- newParams(initConds, startConds, epidemicType)
			if (epidemicType == 3) {
				# Fit k+1 epidemics with new SIR sub epidemic
				eval <- fitInRangeParallel(setSolver(optimMethod, k+1, epiTypesMulti), i, offsetTimes, offsetData, initCondsMulti, initParamsMulti, epiTypesMulti, ts, k+1, c(minTRange:(i-maxTRange)), plotConfig)
			} else if (epidemicType == 1) {
				# Fit k+1 epidemics with set t0 at i
				eval <- fitInRangeParallel(setSolver(optimMethod, k+1, epiTypesMulti), i, offsetTimes, offsetData, initCondsMulti, initParamsMulti, epiTypesMulti, ts, k+1, c((i-1):(i-1)), plotConfig)
			}
			# TODO: Does the residuals array need to be updated here? = Only consider residuals in initial fitInRangeParallel above - will be updated in next loop
			# residuals <- nIncResiduals[(i-(nResiduals-1)):i]
			multiRSquare <- eval$optimRSquare
			# If K + 1 epidemic is better then set k + 1 from now on
			# if (multiRSquare > rSquare) {
				# Set k+1 epidemics from now on
				k <- k + 1
				# Update parameters to continue fitting with k+1 epidemics
				multiParams <- eval$multiParams
				maxt <- eval$optimTime
				ts <- c(ts,maxt)
				rSquare <- multiRSquare
				initConds <- initCondsMulti
				epiTypes <- epiTypesMulti
			# }
		}
		totalRSqaure <- totalRSqaure + rSquare
		# Update the initial parameters for the next fitting
		initParams <- multiParams
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
		meanDiffRes <- myMeanDiff(absResiduals)
		sdDiffRes <- mySdDiff(absResiduals, meanDiffRes)
		diffRes <- abs(residuals[resLength - nRes + 1] - residuals[resLength - nRes])
		print(paste("RSq", rSquare))
		print(paste("meanRes", meanRes))
		print(paste("sdRes ", sdRes))
		# print(paste("MeanDiffRes", meanDiffRes))
		# print(paste("SdDiffRes", sdDiffRes))
		# print(paste("DiffRes", diffRes))
		# If current residual sd is above zero check if last n residuals are above set number of sd
		if (sdRes > 0) {
			minIncRes1 <- Inf
			# Index of first residual to check
			startResIndex <- resLength - nRes + 1
			# Assume incRes is True and check condition for all n residuals
			# sameSign <- max(residuals[(startResIndex + 1):resLength]) < 0 || min(residuals[(startResIndex + 1):resLength]) > 0
			# sigIncRes <- diffRes > (meanDiffRes + 2 * sdDiffRes)
			minIncRes <- min(abs(residuals[startResIndex:resLength]))
			print(minIncRes)
		}
		# print(paste("SameSign",sameSign))
		# Set epidemic type according to residual limit
		# sirSD <- meanRes + (sdRes * 2)
		# spikeSD <- meanRes + (sdRes * 4)
		sirSD <- meanRes + sdRes * 3
		spikeSD <- meanRes + sdRes * 5
		# If minimum residual increase is more than required, then set type
		if (minIncRes > spikeSD) { #&& sameSign) {
			print("Spike Set")
			type <- 1
		} else if (minIncRes > sirSD) { #&& sameSign) {
			print("SIR Set")
			type <- 3
		}
	}
	type
}

newParams <- function(initVec, startVec, epidemicType) {
	if (epidemicType == 1) {
		# Update params for spike epidemic take I0 and gamma from initial
		newVec <- c(initVec, startVec[2])
	} else if (epidemicType == 3) {
		# Update params for SIR epidemic
		newVec <- c(initVec, startVec)
	}

}