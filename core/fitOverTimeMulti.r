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
	# Step size for iterative fitting
	step <- 1
	# Initial t0 value
	ts <- c(1)
	# Set the number of epidemics
	k <- 1

	# All evaluation vector
	evalList <- c()
	
	# Number of increasing residuals
	nResiduals <- 3
	nIncResiduals <- c()
	
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
		nIncResiduals <- c(nIncResiduals, eval$finalResidual)
		# lastNResiduals <- nIncResiduals[(i-(nResiduals-1)):i]
		lastNResiduals <- nIncResiduals


		# Try to improve fit if rSquare has deteriorated
		lim <- thresholds$lim
		diff <- thresholds$diff
		incRes <- incResiduals(lastNResiduals, nResiduals, diff)
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
			# lastNResiduals <- nIncResiduals[(i-(nResiduals-1)):i]
			multiRSquare <- eval$optimRSquare
			# If k+1 is significantly better then continue with k+1 fit
			# if ((multiRSquare - rSquare) > diff || incResiduals(lastNResiduals, nResiduals)) {
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
		# Update the initial parameters for the next fitting
		initParams <- multiParams
	}
	
	# Save all params
	save(evalList, file=plotConfig$dataFile)
	# save.image(plotConfig$envFile)
}

incResiduals <- function(lastNResiduals, n, diff) {
	# # Get the last n residuals
	# lastNResiduals <- (lastNResiduals[!is.na(lastNResiduals)])
	
	# # Initialise check variables
	incResiduals <- FALSE
	# inc <- TRUE
	# sameSign <- TRUE
	# aboveLimit <- (abs(lastNResiduals[1]) > diff)
	# print(abs(lastNResiduals[1])) 

	# if (length(lastNResiduals) == n) {
	# 	for (i in 2:n) {
	# 		# Check magnitude is continuously increasing
	# 		inc <- inc && (abs(lastNResiduals[i]) > abs(lastNResiduals[i-1]))
	# 		# Check residuals are the same sign
	# 		sameSign <- sameSign && ((lastNResiduals[i]<0) == (lastNResiduals[i-1]<0))
	# 		# Check magnitude of residual
	# 		print(abs(lastNResiduals[i])) 
	# 		aboveLimit <- aboveLimit && (abs(lastNResiduals[i]) > diff)
	# 	}
	# 	# Conjunction of check variables indicates new epidemic
	# 	incResiduals <- aboveLimit
	# 	print("IM") 
	# 	print(incResiduals)
	# }
	# incResiduals

	# Calculate mean and standard deviation of all previous residuals
	if (length(lastNResiduals) > 1) {
		# Mean and SD of previous residuals
		# print(lastNResiduals)
		meanRes <- myMean(lastNResiduals[(1:length(lastNResiduals) - 1)])
		sdRes <- mySd(lastNResiduals[1:(length(lastNResiduals) - 1)], meanRes)
		# print("SD: "); print(sdRes)
		# print("MEAN: "); print(meanRes)
		# print("Last Res: "); print(lastNResiduals[length(lastNResiduals)])
		if (abs(lastNResiduals[length(lastNResiduals)]) > 2*sdRes) {
			incResiduals <- TRUE
			print("2SD")
		}
	}
	incResiduals
}