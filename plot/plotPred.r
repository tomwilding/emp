plotPred <- function(times, data, offsets, thresholds, initParams, initConds, plotConfig) {
	# Unpack settings
	minTruncation <- offsets$minTruncation
	startOffset <- offsets$startOffset
	endOffset <- offsets$endOffset
	lim <- thresholds$lim
	evalPreds <- c()
	evalPredsAR <- c()

	# Colours for plotting
	cl <- c("red","cyan","forestgreen","goldenrod2","red4")

	# Loop through all objects

	end <- length(evalList)
	for(i in (minTruncation + 1):end) {
		
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
		evalPrev <- evalList[[i - 1]]
		allEval <- evalPrev$allEval$multiInf

		# Update previous prediction using AR model
		# Get past residuals
		res <- evalPrev$residuals
		# Fit AR model to all past residuals
		arModel <- ar(res)
		nextIncRes <- predict(arModel, n.ahead=1)$pred
		# print(nextIncRes)

		# Plot prediction without AR
		evalPreds[i] <- allEval[i]

		# Plot prediction with AR
		evalPredsAR[i] <- allEval[i] + nextIncRes
	}

	# Calculate SSE
	# SSE of epi
	print(evalPreds[(minTruncation + 1):length(evalPreds)])
	print(evalPredsAR[(minTruncation + 1):length(evalPredsAR)])
	sseEpi <- ssError(evalPreds[(minTruncation + 1):length(evalPreds)], offsetData[(minTruncation + 1):length(offsetData)])
	
	# SSE of AR
	sseAR <- ssError(evalPredsAR[(minTruncation + 1):length(evalPredsAR)], offsetData[(minTruncation + 1):length(offsetData)])

	# SSE of offset data
	# Shift data
	shiftOffsetData <- c(0,offsetData)
	print(shiftOffsetData[(minTruncation + 1):length(offsetData)])
	print(offsetData[(minTruncation + 1):length(offsetData)])
	sseDiff <- ssError(shiftOffsetData[(minTruncation + 1):length(offsetData)], offsetData[(minTruncation + 1):length(offsetData)])

	print(paste("EpiRS", sseEpi))
	print(paste("EpiARRS", sseAR))
	print(paste("ShiftRS", sseDiff))

	# Set graph settings
	setEPS()
	postscript(paste(plotConfig$fileName, "allBlurPrediction", sep=''))	

	# Main plot
	par(mar=c(6.1,4.1,4.1,2.1))
	plot(offsetTimes, offsetData, xlab='Epochs', ylab='Infected Individuals', col='steelblue', type="l")
	title(main=plotConfig$title, cex.main=0.9, cex.axis=0.8)
	daysText <- paste("Epochs after outbreak = ", i)
	mtext(daysText, 3, cex=0.8)

	# Plot actual data point at this time
	lines(offsetTimes, evalPreds, col='red')
	dev.off()
}