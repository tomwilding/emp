plotPred <- function(times, data, offsets, thresholds, initParams, initConds, plotConfig) {
	# Unpack settings
	minTruncation <- offsets$minTruncation
	startOffset <- offsets$startOffset
	endOffset <- offsets$endOffset
	lim <- thresholds$lim
	evalPreds <- c()
	evalPredsAR <- c()
	meanPredsAR <- c()
	futurePredsRS <- c()

	allSSEAR <- c()
	allSSEDiff <- c()
	for (range in 1:20) {
		# Loop through all objects
		predOffset <- 1
		step <- 1
		end <- length(evalList)

		# Take data set within specified offset
		offsetTimes <- times[startOffset:(length(times)-endOffset)]
		offsetData <- data[startOffset:(length(data)-endOffset)]

		# Mean of data
		# meanPred <- myMean(inRangeData)

		# AR model order
		# AROrder <- 4
		quarter <- round(end/4)*2
		# quartiles <- c(quarter, round(quarter*2), round(quarter*3))

		inRangeTimes <- offsetTimes[quarter:(quarter+range)]
		inRangeData <- offsetData[quarter:(quarter+range)]
		# print(quartiles)
		evalPreds <- c()
		evalPredsAR <- c()
		p <- 1
		for(i in seq(from=quarter, to=quarter + range, by=step)) {
			# Plot predicted data point for this time at previous fitting
			evalPrev <- evalList[[quarter - 1]]
			allEval <- evalPrev$allEval$multiInf
			# Update previous prediction using AR model
			# Get past residuals
			prevResiduals <- evalPrev$residuals
			nextIncRes <- 0
			# Fit AR model to all past residuals
			arModel <- ar(prevResiduals)
			# arimaModel <- arima(prevResiduals, order=c(4,0,0))
			nextIncRes <- predict(arModel, n.ahead=range + 1)$pred
			# nextIncRes1 <- predict(arModel, n.ahead=predOffset + 1)$pred	
			# Plot prediction without AR
			# evalPreds[i - 1] <- allEval[(i - 1)]
			# evalPredsAR[i - 1] <- allEval[(i - 1)] + nextIncRes[length(nextIncRes)]
			# meanPredsAR[i - 1] <- meanPred + nextIncRes[length(nextIncRes)]
			evalPreds <- c(evalPreds, allEval[i])
			evalPredsAR <- c(evalPredsAR, allEval[i] + nextIncRes[p])
			# meanPredsAR[i] <- meanPred + nextIncRes1[length(nextIncRes1)] 
			p <- p + 1
		}

		# Calculate SSE
		# SSE of epi
		inRangeEvalPreds <- evalPreds[(minTruncation + predOffset):length(evalPreds)]
		sseEpi <- ssError(inRangeEvalPreds, inRangeData)
		rSqEpi <- rSquareError(inRangeEvalPreds, inRangeData)

		# SSE of AR
		# inRangeEvalPredsAR <- evalPredsAR[(minTruncation + predOffset):length(evalPredsAR)]
		# inRangeEvalPredsAR[101] <- inRangeData[100]
		# inRangeEvalPredsAR[222] <- inRangeData[221]
		# inRangeEvalPredsAR[227] <- inRangeData[226]
		sseAR <- ssError(evalPredsAR, inRangeData)
		rSqAR <- rSquareError(evalPredsAR, inRangeData)
		print(evalPredsAR)
		# SSE of offset data
		# Shift data
		inRangeShiftOffset <- rep(offsetData[(quarter - 1)], range)
		sseDiff <- ssError(inRangeShiftOffset, inRangeData)
		rSqDiff <- rSquareError(inRangeShiftOffset, inRangeData)

		allSSEAR <- c(allSSEAR, sseAR)
		allSSEDiff <- c(allSSEDiff, sseDiff)
	}
	print(allSSEAR)

	plot(c(1:20), allSSEAR, type="l", ylim=c(0, 1e9))
	lines(allSSEDiff, col="steelblue")
	# SSE of mean AR
	# inRangeMeanPredsAR <- meanPredsAR[(minTruncation + predOffset):length(meanPredsAR)]
	# sseMeanAR <- ssError(inRangeMeanPredsAR, inRangeData)
	# rSqMeanAR <- rSquareError(inRangeMeanPredsAR, inRangeData)

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
	# n <- length(inRangeData)
	# sigSqAR <- sseAR / n
	# sigSqDiff <- sseDiff / n
	# aicAR <- log(sigSqAR) + (n + 2*AROrder) / n
	# aicDiff <- log(sigSqDiff) + (n + 2*n) / n

	# n <- length(inRangeDataWR)
	# sigSqARWR <- sseARWR / n
	# sigSqDiffWR <- sseDiffWR / n
	# aicARWR <- log(sigSqARWR) + (n + 2*AROrder) / n
	# aicDiffWR <- log(sigSqDiffWR) + (n + 2*n) / n

	# Median Absolute Difference
	# medEpiDev <- median(abs(inRangeData - inRangeEvalPredsAR))
	# medShiftDev <- median(abs(inRangeData - inRangeShiftOffset))
	# medADAR <- median(abs((inRangeEvalPredsAR - medEpiDev)))
	# medADDiff <- median(abs(inRangeShiftOffset - medShiftDev))

	# madAR <- mad(inRangeEvalPredsAR, inRangeData)
	# madDiff <- mad(inRangeShiftOffset, inRangeData)

	# print(paste("EpiSS", sseEpi))
	# print(paste("EpiARSS", sseAR), quote=FALSE)
	# print(paste("ShiftSS", sseDiff), quote=FALSE)
	# # print(paste("MeanARSS", sseMeanAR), quote=FALSE)

	# print(paste("EpiRS", rSqEpi), quote=FALSE)
	# print(paste("EpiARRS", rSqAR), quote=FALSE)
	# print(paste("ShiftRS", rSqDiff), quote=FALSE)
	# # print(paste("MeanARRS", rSqMeanAR), quote=FALSE)

	# # print(paste("EpiARAIC", aicAR), quote=FALSE)
	# # print(paste("ShiftAIC", aicDiff), quote=FALSE)

	# print(paste("EpiARMedAD", medADAR), quote=FALSE)
	# print(paste("ShiftMedAD", medADDiff), quote=FALSE)

	# print(paste("EpiARMAD", madAR), quote=FALSE)
	# print(paste("ShiftMAD", madDiff), quote=FALSE)	

	# # print(paste("EpiMean", myMean(abs(inRangeEvalPreds - inRangeData))))
	# # print(paste("EpiARMean", myMean(abs(inRangeEvalPredsAR - inRangeData))))	
	# # print(paste("ShiftMean", myMean(abs(inRangeShiftOffset - inRangeData))))

	# # Set graph settings
	# setEPS()
	# postscript(paste(plotConfig$fileName, "allCallPredictionAR.eps", sep=''))	

	# # Main plot
	# par(mar=c(6.1,4.1,4.1,2.1))
	# plot(inRangeTimes, inRangeData, xlab='Time (Days)', ylab='Predicted Infected Individuals', col='steelblue', type="l")
	# title(main="Next Day Predictions of Synthedemic model with AR Residual Refinement", cex.main=0.9, cex.axis=0.8)
	# # daysText <- paste("Day", i)
	# # mtext(daysText, 3, cex=0.8)

	# # Plot actual data point at this time
	# print(length(inRangeTimes))
	# print(length(inRangeEvalPredsAR))
	# print(length(inRangeEvalPreds))
	# lines(inRangeTimes, inRangeEvalPredsAR, col='red')
	# lines(inRangeTimes, inRangeEvalPreds, col='black')
	# # lines(inRangeTimes, inRangeMeanPredsAR, col='blue')
	# # Legend
	# legendText <- c("Data", "Synthedemic", "Synthedemic with AR Refinement")
	# lineType <- c(1, 1, 1)
	# col <- c("steelblue", "black", "red")
	# legend("topright",legendText, col=col, lty=lineType, cex=0.8)
	# dev.off()
}