plotRSq <- function(times, data, offset, thresholds, initParams, initConds, plotConfig) {
	require("epi")

	# Unpack settings
	minTruncation <- offsets$minTruncation
	startOffset <- offsets$startOffset
	endOffset <- offsets$endOffset
	lim <- thresholds$lim
	diff <- thresholds$diff

	# Array to hold all rSqaure evaluations
	allRSqaurePast <- c()
	allRSqaureNF <- c()

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
		# Future RSquare of Near Future to evaluate over
		NF <- 7
		# Calculate RSqaure
		rSquarePast <- eval$optimRSquare
		rSquareNF <- rSquareError(allEval$multiInf[1:i+NF], offsetData[1:i+NF])

		allRSqaurePast <- c(allRSqaurePast, rSquarePast)
		allRSqaureNF <- c(allRSqaureNF, rSquareNF)
	}

	# Plot past rSqaure over time
	# Set graph settings
	setEPS()
	graphName <- "rSquarePast.eps"
	postscript(paste(plotConfig$fileName, graphName, sep=''), height=3)

	plot(c(minTruncation:length(offsetData)), allRSqaurePast, xlab='Time (Days)', ylab='Past RSquare Error', col='steelblue', type="l")
	title("Past RSquare Error at each Day for Blurred Lines")

	dev.off()

}