plotRSq <- function(times, data, offset, thresholds, initParams, initConds, plotConfig) {
	# Unpack settings
	minTruncation <- offsets$minTruncation
	startOffset <- offsets$startOffset
	endOffset <- offsets$endOffset
	lim <- thresholds$lim
	diff <- thresholds$diff

	# Array to hold all rSqaure evaluations
	allRSqaurePast <- c()
	allRSqaureT7 <- c()

	# Loop through all objects
	for(i in minTruncation:length(evalList)) {	
		
		# Take data set within specified offset
		offsetTimes <- times[startOffset:(length(times)-endOffset)]
		offsetData <- data[startOffset:(length(data)-endOffset)]
		# Truncate data set
		truncTimes <- offsetTimes[1:i]
		truncData <- offsetData[1:i]

		# Fine Times for evaluation
		timeStep <- 0.05
		fineTimes <- breakTime(offsetTimes, timeStep)

		# Get graph object for this iteration
		eval <- evalList[[i]]
		allEval <- eval$allEval
		allEvalFine <- eval$allEvalFine

		# RSquare and labels
		# Future RSqaure window to evaluate over
		T7 <- 7
		# Calculate RSqaure
		rSquarePast <- eval$optimRSquare
		rSquareT7 <- rSquareError(allEval$multiInf[1:i+T7], offsetData[1:i+T7])

		allRSqaurePast <- c(allRSqaurePast, rSquarePast)
		allRSqaureT7 <- c(allRSqaureT7, rSquareT7)
	}

	# Plot past rSqaure over time
	# Set graph settings
	setEPS()
	graphName <- "rSqaurePast.eps"
	postscript(paste(plotConfig$fileName, graphName, sep=''), height=3)

	plot(c(1:(length(offsetData)+1 - minTruncation)), allRSqaurePast, xlab='Epochs', ylab='Past RSqaure Error', col='steelblue')
	title("Past RSquare Error at each Epoch")
	# Past RSq
	points(c(1:(length(offsetData)+1 - minTruncation)), allRSqaurePast, col='steelblue')
	lines(c(1:(length(offsetData)+1 - minTruncation)), allRSqaurePast, col="steelblue")

	dev.off()

}