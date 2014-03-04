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
	nResiduals <- 2
	residuals <- c()
	
	################################################# Decompose Epidemics ################################################
	# Truncate the data to i data points from minTruncation within offset data
	for (i in seq(from=minTruncation, to=maxTruncation, by=step)) {
		# Fit k epidemics
		print("### Fit k", quote=FALSE); print(paste(c("fitting "," of "), c(i, maxTruncation)), quote=FALSE)
		eval <- fitInRangeParallel(setSolver(optimMethod, k, epiTypes), i, offsetTimes, offsetData, initConds, initParams, epiTypes, ts, k, c(minTRange:(i-maxTRange)), plotConfig)
		# Store eval
		evalList[[i]] <- eval

		maxt <- eval$optimTime
		ts[k] <- maxt
		rSquare <- eval$optimRSquare
		multiParams <- eval$multiParams
		# Get last residual and update last n residuals vector
		residuals <- c(residuals, eval$finalResidual)


		# Try to improve fit if rSquare has deteriorated
		lim <- thresholds$lim
		diff <- thresholds$diff
		epidemicType <- getEpidemicType(residuals, nResiduals)
		if ((rSquare < lim) && (epidemicType > 0)) {
			# Set new epidemic type in epidemic type array
			epiTypes <- c(epiTypes, epidemicType)
			# Try k+1 epidemics
			print(">>> Fit k+1", quote=FALSE)
			# Set initial parameters
			initParamsMulti <- newParams(initParams, startParams, epidemicType)
			initCondsMulti <- newParams(initConds, startConds, epidemicType)
			# Fit k+1 epidemics
			eval <- fitInRangeParallel(setSolver(optimMethod, k+1, epiTypes), i, offsetTimes, offsetData, initCondsMulti, initParamsMulti, epiTypes, ts, k+1, c(minTRange:(i-maxTRange)), plotConfig)
			# TODO: Does the residuals array need to be updated here? = Only consider residuals in initial fitInRangeParallel above - will be updated in next loop
			# residuals <- nIncResiduals[(i-(nResiduals-1)):i]
			multiRSquare <- eval$optimRSquare
			# If k+1 is significantly better then continue with k+1 fit
			print("!!!Set k+1", quote=FALSE)
			# Set k+1 epidemics from now on
			k <- k + 1
			# Update parameters to continue fitting with k+1 epidemics
			multiParams <- eval$multiParams
			maxt <- eval$optimTime
			ts <- c(ts,maxt)
			rSquare <- multiRSquare
			initConds <- initCondsMulti
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
getEpidemicType <- function(residuals, nRes) {
	incRes <- c()
	type <- 0
	# Standard deviation of residuals
	sdRes <- 0
	# Get the last n residuals
	residuals <- (residuals[!is.na(residuals)])
	# If less then one residual then set residual check to false
	minIncRes <- Inf
	# Ensure more than one residual before the last n residuals to calculate sdRes
	if (length(residuals) > nRes + 1) {
		# Get standard deviation of residuals before the ones considered
		meanRes <- myMean(residuals[(1:(length(residuals) - nRes))])
		sdRes <- mySd(residuals[1:(length(residuals) - nRes)], meanRes)
		print(paste("sdRes ", sdRes))
		# If current residual sd is above zero check if last n residuals are above set number of sd
		if (sdRes > 0) {
			# Assume incRes is True and check condition for all n residuals
			incResiduals <- TRUE
			for (i in 1:nRes) {
				# Check magnitude of residuals
				startResIndex <- length(residuals) - nRes
				print(residuals[startResIndex + i])
				minIncRes <- min(minIncRes, (abs(residuals[startResIndex + i])))
			}
		}
		# Set epidemic type according to residual limit
		sirSD <- sdRes
		spikeSD <- sdRes * 8
		# If minimum residual increase is more than required, then set type
		if (minIncRes > spikeSD) {
			type <- 1
		}
		else if (minIncRes > sirSD * 2) {
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