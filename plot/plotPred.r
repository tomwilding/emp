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

	# Take data set within specified offset
	offsetTimes <- times[startOffset:(length(times)-endOffset)]
	offsetData <- data[startOffset:(length(data)-endOffset)]

	# Loop through all objects
	end <- length(evalList)
	predOffset <- 1
	for(i in (minTruncation + predOffset):end) {

		# Plot predicted data point for this time at previous fitting
		evalPrev <- evalList[[i - predOffset]]
		allEval <- evalPrev$allEval$multiInf
		
		# Update previous prediction using AR model
		# Get past residuals
		# lastEpiTime <- evalPrev$optimTimes[length(evalPrev$optimTimes)]
		# prevResiduals <- evalPrev$residuals[lastEpiTime:i-1]
		prevResiduals <- evalPrev$residuals
		# Fit AR model to all past residuals
		arModel <- ar(prevResiduals)
		nextIncRes <- predict(arModel, n.ahead=predOffset)$pred

		# Plot prediction without AR
		evalPreds[i] <- allEval[i]

		# Plot prediction with AR
		evalPredsAR[i] <- allEval[i] + nextIncRes[predOffset]
	}

	# Calculate SSE
	# SSE of epi
	inRangeEvalPreds <- evalPreds[(minTruncation + predOffset):length(evalPreds)]
	inRangeData <- offsetData[(minTruncation + predOffset):length(offsetData)]
	inRangeTimes <- offsetTimes[(minTruncation + predOffset):length(offsetTimes)]
	sseEpi <- ssError(inRangeEvalPreds, inRangeData)
	rSqEpi <- rSquareError(inRangeEvalPreds, inRangeData)

	# SSE of AR
	inRangeEvalPredsAR <- evalPredsAR[(minTruncation + predOffset):length(evalPredsAR)]
	sseAR <- ssError(inRangeEvalPredsAR, inRangeData)
	rSqAR <- rSquareError(inRangeEvalPredsAR, inRangeData)

	# print(inRangeEvalPredsAR - inRangeData)
	# print(inRangeData)

	# SSE of offset data
	# Shift data
	offsetPadding <- numeric(predOffset)
	shiftOffsetData <- c(offsetPadding,inRangeData)
	inRangeShiftOffset <- shiftOffsetData[1:length(inRangeData)]
	sseDiff <- ssError(inRangeShiftOffset, inRangeData)
	rSqDiff <- rSquareError(inRangeShiftOffset, inRangeData)

	# Without repeat data points
	inRangeShiftOffsetWithoutRepeats <- c(shiftOffsetData[1:115], shiftOffsetData[122:length(inRangeData)])
	inRangeDataWithoutRepeats <- c(inRangeData[1:114], inRangeData[121:length(inRangeData)])
	inRangeEvalPredsARWR <- c(inRangeEvalPredsAR[1:114], inRangeEvalPredsAR[121:length(inRangeEvalPredsAR)])
	sseDiffWR <- ssError(inRangeShiftOffsetWithoutRepeats, inRangeDataWithoutRepeats)
	sseARWR <- ssError(inRangeEvalPredsARWR, inRangeDataWithoutRepeats)
	rSqARWR <- rSquareError(inRangeEvalPredsARWR, inRangeDataWithoutRepeats)
	rSqDiffWR <- rSquareError(inRangeShiftOffsetWithoutRepeats, inRangeDataWithoutRepeats)
	# print(inRangeEvalPredsARWR - inRangeDataWithoutRepeats)
	# print(inRangeShiftOffsetWithoutRepeats - inRangeDataWithoutRepeats)

	# print(paste("EpiSS", sseEpi))
	print(paste("EpiARSSWR", sseARWR))
	print(paste("ShiftSSWR", sseDiffWR))
	# print(paste("EpiRS", rSqEpi))
	print(paste("EpiARRSWR", rSqARWR))
	print(paste("ShiftRSWR", rSqDiffWR))
	# print(paste("EpiMean", myMean(abs(inRangeEvalPreds - inRangeData))))
	# print(paste("EpiARMean", myMean(abs(inRangeEvalPredsAR - inRangeData))))	
	# print(paste("ShiftMean", myMean(abs(inRangeShiftOffset - inRangeData))))

	# Set graph settings
	setEPS()
	postscript(paste(plotConfig$fileName, "allBlurPredictionAR", sep=''))	

	# Main plot
	par(mar=c(6.1,4.1,4.1,2.1))
	plot(inRangeTimes, inRangeData, xlab='Epochs', ylab='Infected Individuals', col='steelblue', type="l")
	title(main=plotConfig$title, cex.main=0.9, cex.axis=0.8)
	daysText <- paste("Epochs after outbreak = ", i)
	mtext(daysText, 3, cex=0.8)

	# Plot actual data point at this time
	lines(inRangeTimes, inRangeEvalPredsAR, col='red')
	abline(v=106)
	dev.off()
}