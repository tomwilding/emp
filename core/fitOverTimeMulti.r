fitOverTimeMulti <- function(optimMethod, times, data, initConds, initParams, offsets, thresholds, plotConfig) {
	
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
		eval <- fitInRangeParallel(setSolver(optimMethod, k), i, offsetTimes, offsetData, initConds, initParams, ts, k, c(minTRange:(i-maxTRange)), 1, plotConfig)
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
		incRes <- incResiduals(residuals, nResiduals)
		if ((rSquare < lim) && incRes) {
		# if (rSquare < lim) {
			# Try k+1 epidemics
			print(">>> Fit k+1", quote=FALSE)
			# Set initial parameters
			initParamsMulti <- c(initParams, startParams)
			initCondsMulti <- c(initConds, startConds)
			# Fit k+1 epidemics
			eval <- fitInRangeParallel(setSolver(optimMethod, k+1), i, offsetTimes, offsetData, initCondsMulti, initParamsMulti, ts, k+1, c(minTRange:(i-maxTRange)), 2, plotConfig)
			# TODO: Does the residuals array need to be updated here?
			# residuals <- nIncResiduals[(i-(nResiduals-1)):i]
			multiRSquare <- eval$optimRSquare
			# If k+1 is significantly better then continue with k+1 fit
			# if ((multiRSquare - rSquare) > diff || incResiduals(residuals, nResiduals)) {
			# if ((multiRSquare - rSquare) > diff) {
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

incResiduals <- function(residuals, n) {
	# # Get the last n residuals
	residuals <- (residuals[!is.na(residuals)])
	
	# # Initialise check variables
	incResiduals <- TRUE
	print(length(residuals)); print(n)
	if (length(residuals) > 1) {
		meanRes <- myMean(residuals[(1:length(residuals) - 1)])
		sdRes <- mySd(residuals[1:(length(residuals) - 1)], meanRes)
		print(sdRes)
		for (i in 1:n) {
			# Check magnitude of residuals
			minResIndex <- length(residuals) - n
			print(residuals[minResIndex + i])
			incResiduals <- incResiduals && (abs(residuals[minResIndex + i]) > sdRes)
			# readline()
		}
		print(incResiduals)
	}
	incResiduals <- FALSE
}