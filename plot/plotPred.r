plotPred <- function(times, data, offsets, thresholds, initParams, initConds, plotConfig) {
	# Unpack settings
	minTruncation <- offsets$minTruncation
	startOffset <- offsets$startOffset
	endOffset <- offsets$endOffset
	lim <- thresholds$lim
	evalPreds <- c()
	evalPredsAR <- c()
	meanPredsAR <- c()

	# Loop through all objects
	predOffset <- 1
	end <- length(evalList)

	# Take data set within specified offset
	offsetTimes <- times[startOffset:(length(times)-endOffset)]
	offsetData <- data[startOffset:(length(data)-endOffset)]
	inRangeTimes <- offsetTimes[(minTruncation + predOffset):length(offsetTimes)]
	inRangeData <- offsetData[(minTruncation + predOffset):end]

	# Mean of data
	meanPred <- myMean(inRangeData)

	# AR model order
	AROrder <- 4

	# Residual window
	# window <- 4

	for(i in (minTruncation + predOffset):end) {
		# Plot predicted data point for this time at previous fitting
		evalPrev <- evalList[[i - predOffset]]
		allEval <- evalPrev$allEval$multiInf
		
		# Update previous prediction using AR model
		# Get past residuals
		# lastEpiTime <- evalPrev$optimTimes[length(evalPrev$optimTimes)]
		# prevResiduals <- evalPrev$residuals[lastEpiTime:length(evalPrev$residuals)]
		prevResiduals <- evalPrev$residuals
		# prevResidualsWindow <- prevResiduals[length(prevResiduals)-window:length(prevResiduals)]
		# Fit AR model to all past residuals
		arModel <- ar(prevResiduals, FALSE, 4)
		# arimaModel <- arima(prevResiduals, order=c(4,0,0))
		nextIncRes <- predict(arModel, n.ahead=predOffset)$pred
		# Plot prediction without AR
		evalPreds[i] <- allEval[i]

		# Plot prediction with AR
		# for (j in 0:(length(nextIncRes) - 1)) {
		# 	evalPredsAR[i - j] <- allEval[i - j] + nextIncRes[length(nextIncRes) - j]
		# }
		evalPredsAR[i] <- allEval[i] + nextIncRes[length(nextIncRes)]
		meanPredsAR[i] <- meanPred + nextIncRes[length(nextIncRes)]
	}
	# print(arModel)
	# print(arimaModel)

	# Calculate SSE
	# SSE of epi
	inRangeEvalPreds <- evalPreds[(minTruncation + predOffset):length(evalPreds)]
	sseEpi <- ssError(inRangeEvalPreds, inRangeData)
	rSqEpi <- rSquareError(inRangeEvalPreds, inRangeData)

	# SSE of AR
	inRangeEvalPredsAR <- evalPredsAR[(minTruncation + predOffset):length(evalPredsAR)]
	# inRangeEvalPredsAR[101] <- inRangeData[100]
	# inRangeEvalPredsAR[222] <- inRangeData[221]
	# inRangeEvalPredsAR[227] <- inRangeData[226]
	sseAR <- ssError(inRangeEvalPredsAR, inRangeData)
	rSqAR <- rSquareError(inRangeEvalPredsAR, inRangeData)

	# SSE of offset data
	# Shift data
	offsetPadding <- numeric(predOffset)
	shiftOffsetData <- c(offsetPadding,inRangeData)
	inRangeShiftOffset <- shiftOffsetData[1:length(inRangeData)]
	sseDiff <- ssError(inRangeShiftOffset, inRangeData)
	rSqDiff <- rSquareError(inRangeShiftOffset, inRangeData)

	# SSE of mean AR
	inRangeMeanPredsAR <- meanPredsAR[(minTruncation + predOffset):length(meanPredsAR)]
	sseMeanAR <- ssError(inRangeMeanPredsAR, inRangeData)
	rSqMeanAR <- rSquareError(inRangeMeanPredsAR, inRangeData)

	# Without repeat data points
	# inRangeShiftOffsetWR <- c(shiftOffsetData[1:115], shiftOffsetData[122:length(inRangeData)])
	# inRangeDataWR <- c(inRangeData[1:(115-predOffset)], inRangeData[(122-predOffset):length(inRangeData)])
	# inRangeEvalPredsWR <- c(inRangeEvalPreds[1:(115-predOffset)], inRangeEvalPreds[(122-predOffset):length(inRangeEvalPreds)])
	# inRangeEvalPredsARWR <- c(inRangeEvalPredsAR[1:(115-predOffset)], inRangeEvalPredsAR[(122-predOffset):length(inRangeEvalPredsAR)])
	# sseEpiWR <- ssError(inRangeEvalPredsWR, inRangeDataWR)
	# sseARWR <- ssError(inRangeEvalPredsARWR, inRangeDataWR)
	# sseDiffWR <- ssError(inRangeShiftOffsetWR, inRangeDataWR)
	# rSqARWR <- rSquareError(inRangeEvalPredsARWR, inRangeDataWR)
	# rSqDiffWR <- rSquareError(inRangeShiftOffsetWR, inRangeDataWR)

	# Calculate Akaike
	n <- length(inRangeData)
	sigSqAR <- sseAR / n
	sigSqDiff <- sseDiff / n
	aicAR <- log(sigSqAR) + (n + 2*AROrder) / n
	aicDiff <- log(sigSqDiff) + (n + 2*n) / n

	# n <- length(inRangeDataWR)
	# sigSqARWR <- sseARWR / n
	# sigSqDiffWR <- sseDiffWR / n
	# aicARWR <- log(sigSqARWR) + (n + 2*AROrder) / n
	# aicDiffWR <- log(sigSqDiffWR) + (n + 2*n) / n

	print(paste("EpiSS", sseEpi))
	print(paste("EpiARSS", sseAR), quote=FALSE)
	print(paste("ShiftSS", sseDiff), quote=FALSE)
	print(paste("MeanARSS", sseMeanAR), quote=FALSE)

	print(paste("EpiRS", rSqEpi), quote=FALSE)
	print(paste("EpiARRS", rSqAR), quote=FALSE)
	print(paste("ShiftRS", rSqDiff), quote=FALSE)
	print(paste("MeanARRS", rSqMeanAR), quote=FALSE)
	
	print(paste("EpiARAIC", aicAR), quote=FALSE)
	print(paste("ShiftAIC", aicDiff), quote=FALSE)


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
	# lines(inRangeTimes, inRangeEvalPreds, col='green')
	# lines(inRangeTimes, inRangeMeanPredsAR, col='blue')
	# abline(v=106)
	dev.off()
}