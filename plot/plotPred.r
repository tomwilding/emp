plotPred <- function(times, data, offsets, thresholds, initParams, initConds, plotConfig) {
	# Unpack settings
	minTruncation <- offsets$minTruncation
	startOffset <- offsets$startOffset
	endOffset <- offsets$endOffset
	lim <- thresholds$lim
	evalPreds <- c()

	# Colours for plotting
	cl <- c("red","cyan","forestgreen","goldenrod2","red4")

	# Loop through all objects

	end <- length(evalList)
	for(i in minTruncation:end) {
		
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
		# eval <- evalList[[i]]
		# allEval <- eval$allEval
		# allEvalFine <- eval$allEvalFine

		# Plot predicted data point for this time at previous fitting
		evalPrev <- evalList[[i]]
		allEval <- evalPrev$allEval$multiInf

		evalPreds[i] <- allEval[i]
	}

	# # Set graph settings
	# setEPS()
	# graphName <- paste("t", i, sep='')
	# graphName <- paste(graphName, ".eps", sep='')
	# postscript(paste(plotConfig$fileName, graphName, sep=''))	

	# Main plot
	par(mar=c(6.1,4.1,4.1,2.1))
	plot(offsetTimes, offsetData, xlab='Epochs', ylab='Infected Individuals', col='steelblue', type="l")
	title(main=plotConfig$title, cex.main=0.9, cex.axis=0.8)
	daysText <- paste("Epochs after outbreak = ", i)
	mtext(daysText, 3, cex=0.8)

	# Plot actual data point at this time
	lines(offsetTimes, evalPreds, col='red')
}