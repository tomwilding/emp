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

	# Step size for iterative fitting
	step <- 1
	# Initial t0 value
	ts <- c(1)

	# All evaluation vector
	evalList <- c()
	
	################################################# Decompose Epidemics ################################################
	# Truncate the data to i data points from 20 within offset data
	for (i in seq(from=minTruncation, to=maxTruncation, by=step)) {
		# Fit k epidemics
		print("------------------------------------------------", quote=FALSE)
		print(paste(c("fitting "," of "), c(i, maxTruncation)), quote=FALSE);
		# Determine epidemic type and fit over required range
		eval <- fitInRangeParallel(setSolver(optimMethod, 1, epiTypes), i, offsetTimes, offsetData, initConds, initParams, epiTypes, ts, 1, c(ts[1]:ts[1]), plotConfig, 1)
		# Update parameters
		maxt <- eval$optimTime
		ts[1] <- maxt
		rSquare <- eval$optimRSquare
		optimParams <- eval$multiParams
		optimConds <- eval$initConds
		print("optimParams")
		print(optimParams)
		print(paste("rs", rSquare))
		# Get last residual and update residuals vector
		residuals <- eval$residuals
		evalList[[i]] <- eval
	}
	
	# Save all params
	save(evalList, file=plotConfig$dataFile)
	# save.image(plotConfig$envFile)
}
