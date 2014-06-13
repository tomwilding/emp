plotPred <- function(times, data, offsets, thresholds, initParams, initConds, plotConfig) {
	require("forecast")
	# Unpack settings
	minTruncation <- offsets$minTruncation
	startOffset <- offsets$startOffset
	endOffset <- offsets$endOffset
	lim <- thresholds$lim
	evalPreds <- c()
	evalPredsAR <- c()
	meanPredsAR <- c()
	futurePredsRS <- c()

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

	# hw <- HoltWinters(ts(inRangeData, frequency=24))
	# print(hw$fitted[,1])
	# plot(hw)
	# plot(fit(hw))
	# print(length(hw$fitted[,1]))
	# print(length(inRangeData))
	# print(rSquareError(c(numeric(24),hw$fitted[,1]), inRangeData))
	# readline()
	allRS <- c()
	allSSE <- c()
	allRSDiff <- c()
	allSSEDiff <- c()
	for (predOffset in c(1:4)) {
		inRangeTimes <- offsetTimes[(minTruncation + predOffset):end]
		inRangeData <- offsetData[(minTruncation + predOffset):end]
		evalPreds <- c()
		evalPredsAR <- c()
		evalPredsARIMA <- c()
		inRangeStepData <- c()
		for(i in seq(from=(minTruncation + (predOffset*step)), to=end, by=step)) {
			# Plot predicted data point for this time at previous fitting
			evalPrev <- evalList[[i - predOffset]]
			allEval <- evalPrev$allEval$multiInf
			# Update previous prediction using AR model
			# Get past residuals
			prevResiduals <- evalPrev$residuals
			# Fit AR model to all past residuals
			nextIncRes <- 0
			nextARIMAPred <- 0
			# order <- 4
			# if (length(prevResiduals) > order + predOffset) {
			arModel <- ar(prevResiduals)
			nextIncRes <- predict(arModel, n.ahead=predOffset)$pred	

				# ARIMA model
				# prevData <- diff(offsetData[1:(i - predOffset)], differences=2)
				# prevData <- offsetData[1:(i - predOffset)]
				# arimaModel <- arima(prevData, c(order,2,1))
				# nextARIMAPred <- predict(arimaModel, n.ahead=step)$pred
			# }
			# evalPreds <- c(evalPreds, allEval[i-1])
			evalPreds <- c(evalPreds, allEval[i])
			# evalPredsAR <- c(evalPredsAR, allEval[i-1] + nextIncRes[length(nextIncRes)-1])
			evalPredsAR <- c(evalPredsAR, allEval[i] + nextIncRes[length(nextIncRes)])
			evalPredsARIMA <- c(evalPredsARIMA, nextARIMAPred)
		}
		allRS <- c(allRS, rSquareError(evalPredsAR, inRangeData))
		allSSE <- c(allSSE, ssError(evalPredsAR, inRangeData))

		# Shift data
		offsetPadding <- numeric(predOffset)
		shiftOffsetData <- c(offsetPadding,inRangeData)
		inRangeShiftOffset <- shiftOffsetData[1:length(inRangeData)]
		allRSDiff <- c(allRSDiff, rSquareError(inRangeShiftOffset, inRangeData))
		allSSEDiff <- c(allSSEDiff, ssError(inRangeShiftOffset, inRangeData))
		
		plot(inRangeTimes, inRangeData, xlab='Time (Days)', ylab='Predicted Infected Individuals', col='steelblue', type="l")
		lines(inRangeTimes, evalPredsAR, col="red")
	}

	print(allRS)
	print(allRSDiff)
	# plot(c(1:4), allRS, type="l")
	# lines(allRSDiff, col="steelblue")
	readline()
	# print(arModel)
	# print(arimaModel)

	# Calculate SSE
	# SSE of epi
	# evalPreds <- evalPreds[(minTruncation + predOffset):length(evalPreds)]
	# sseEpi <- ssError(evalPreds, inRangeData)
	# rSqEpi <- rSquareError(evalPreds, inRangeData)
	# madEpi <- mad(evalPreds, inRangeData)
	# mapeEpi <- mape(evalPreds, inRangeData)
	# rmseEpi <- rmse(evalPreds, inRangeData)
	# raeEpi <- rae(evalPreds, inRangeData)

	# sseAR <- ssError(evalPredsAR, inRangeData)
	# rSqAR <- rSquareError(evalPredsAR, inRangeData)
	# madAR <- mad(evalPredsAR, inRangeData)
	# mapeAR <- mape(evalPredsAR, inRangeData)
	# rmseAR <- rmse(evalPredsAR, inRangeData)
	# raeAR <- rae(evalPredsAR, inRangeData)
	

	# SSE of AR
	# evalPredsAR <- evalPredsAR[(minTruncation + predOffset):length(evalPredsAR)]
	# evalPredsAR[101] <- inRangeData[100]
	# evalPredsAR[222] <- inRangeData[221]
	# evalPredsAR[227] <- inRangeData[226]
	# sseAR <- ssError(evalPredsAR, inRangeData)
	# rSqAR <- rSquareError(evalPredsAR, inRangeData)

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

	# MEDIAN Absolute Difference
	# m <- median(inRangeData)
	# medADAR <- median(abs((evalPredsAR - m)))
	# medADDiff <- median(abs(inRangeShiftOffset - m))


	# rSqARIMA <- rSquareError(evalPredsARIMA, inRangeData)
	# sseARIMA <- ssError(evalPredsARIMA, inRangeData)
	# # Without repeat data points
	# # inRangeShiftOffsetWR <- c(shiftOffsetData[1:141], shiftOffsetData[149:length(inRangeData)])
	# # inRangeDataWR <- c(inRangeData[1:141], inRangeData[149:length(inRangeData)])
	# # evalPredsWR <- c(evalPreds[1:141], evalPreds[149:length(evalPreds)])
	# # evalPredsARWR <- c(evalPredsAR[1:141], evalPredsAR[149:length(evalPredsAR)])
	# # sseEpiWR <- ssError(evalPredsWR, inRangeDataWR)
	# # sseARWR <- ssError(evalPredsARWR, inRangeDataWR)
	# # sseDiffWR <- ssError(inRangeShiftOffsetWR, inRangeDataWR)
	# # rSqARWR <- rSquareError(evalPredsARWR, inRangeDataWR)
	# # rSqDiffWR <- rSquareError(inRangeShiftOffsetWR, inRangeDataWR)

	# print(paste("EpiRS", rSqEpi), quote=FALSE)
	# print(paste("EpiARRS", rSqAR), quote=FALSE)
	# print(paste("ShiftRS", rSqDiff), quote=FALSE)
	# print(paste("ARIMARS", rSqARIMA), quote=FALSE)
	# # print(paste("EpiARRSWR", rSqARWR), quote=FALSE)
	# # print(paste("ShiftRSWR", rSqDiffWR), quote=FALSE)
	# # print(paste("MeanARRS", rSqMeanAR), quote=FALSE)

	# print(paste("EpiSS", sseEpi), quote=FALSE)
	# print(paste("EpiARSS", sseAR), quote=FALSE)
	# print(paste("ShiftSS", sseDiff), quote=FALSE)
	# print(paste("ARIMASS", sseARIMA), quote=FALSE)
	# # print(paste("EpiARSSWR", sseARWR), quote=FALSE)
	# # print(paste("ShiftSSWR", sseDiffWR), quote=FALSE)
	# # print(paste("MeanARSS", sseMeanAR), quote=FALSE)
	
	# # print(paste("EpiARAIC", aicAR), quote=FALSE)
	# # print(paste("ShiftAIC", aicDiff), quote=FALSE)

	# # print(paste("EpiARMedAD", madAR), quote=FALSE)
	# # print(paste("ShiftMedAD", madDiff), quote=FALSE)

	# print(paste("EpiMAD", madEpi), quote=FALSE)
	# print(paste("EpiARMAD", madAR), quote=FALSE)
	# print(paste("ShiftMAD", madDiff), quote=FALSE)

	# print(paste("EpiARMAPE", mapeAR), quote=FALSE)


	# print(paste("EpiARRMSE", rmseAR), quote=FALSE)	
	# # print(paste("EpiARRMSE", rmseAR), quote=FALSE)
	# # print(paste("ShiftRMSE", rmseDiff), quote=FALSE)

	# print(paste("EpiARRAE", raeAR), quote=FALSE)	

	# print(paste("EpiMean", myMean(abs(evalPreds - inRangeData))))
	# print(paste("EpiARMean", myMean(abs(evalPredsAR - inRangeData))))	
	# print(paste("ShiftMean", myMean(abs(inRangeShiftOffset - inRangeData))))

	# # Set graph settings
	# setEPS()
	# postscript(paste(plotConfig$fileName, "allMixPredictionAR.eps", sep=''))	

	# # Main plot
	# par(mar=c(6.1,4.1,4.1,2.1))
	# plot(inRangeTimes, inRangeData, xlab='Time (Days)', ylab='Predicted Infected Individuals', col='steelblue', type="l")
	# title(main="Next Day Predictions of Synthedemic model with AR Residual Refinement", cex.main=0.9, cex.axis=0.8)
	# # daysText <- paste("Day", i)
	# # mtext(daysText, 3, cex=0.8)

	# # Plot actual data point at this time
	# lines(inRangeTimes, evalPredsAR, col='red')
	# lines(inRangeTimes, evalPredsARIMA, col='red')
	# lines(inRangeTimes, inRangeMeanPredsAR, col='blue')
	# # Legend
	# legendText <- c("Data", "Synthedemic", "Synthedemic with AR Refinement")
	# lineType <- c(1, 1, 1)
	# col <- c("steelblue", "black", "red")
	# legend("topright",legendText, col=col, lty=lineType, cex=0.8)
	# dev.off()
}