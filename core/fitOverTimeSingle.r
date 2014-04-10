fitOverTimeSingle <- function(optimMethod, times, data, initConds, initParams, epiTypes, offsets, thresholds, plotConfig) {
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
	window <- 20
	residuals <- c()

	# At the start the start time is set to 0
	startTimePrev <- ts[1]
	startTimeCount <- 10
	repeatStartTimes <- 5

	epiList <- numeric(minTruncation)
	
	################################################# Decompose Epidemics ################################################
	# Truncate the data to i data points from minTruncation within offset data
	for (i in seq(from=minTruncation, to=maxTruncation, by=step)) {
		# Fit k epidemics
		print("### Fit k", quote=FALSE); print(paste(c("fitting "," of "), c(i, maxTruncation)), quote=FALSE)
		# print(paste("Beta", exp(initParams[2])))
		epidemicType <- epiTypes[k]
		# Determine if current epidemic start time is set
		startTimeCount <- countStartTime(ts[k], startTimePrev, startTimeCount)
		startTimePrev <- ts[k]
		print("Count")
		print(startTimeCount)
		# Determine epidemic type and fit over required range
		eval <- fitInRangeParallel(setSolver(optimMethod, k, epiTypes), i, offsetTimes, offsetData, initConds, initParams, epiTypes, ts, k, c(ts[k]:ts[k]), plotConfig)
		maxt <- eval$optimTime
		ts[k] <- maxt
		rSquare <- eval$optimRSquare
		multiParams <- eval$multiParams
		# Get last residual and update residuals vector
		residuals <- c(residuals, eval$finalResidual)

		# TODO: Store eval starts at i offset
		evalList[[i]] <- eval
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