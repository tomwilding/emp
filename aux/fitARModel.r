fitARModel <- function(times, data, offsets, plotConfig) {
	
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

	# All evaluation vector
	evalList <- c()
	
	################################################# AR Iterative Fitting ################################################
	# Truncate the data to i data points from 20 within offset data
	for (i in seq(from=minTruncation, to=maxTruncation, by=step)) {
		truncTimes <- offsetTimes[1:i]
		truncData <- offsetData[1:i]
		arModel <- ar(truncData)
		print(arModel)
		plot(times, data)
		lines(predict(arModel,n.ahead=(maxTruncation - i))$pred,lw=2,col="black")

	}
}