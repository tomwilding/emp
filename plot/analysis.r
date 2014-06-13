analysis <- function(times, data, offsets, thresholds, initParams, initConds, plotConfig) {
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
	step <- 1
	end <- length(evalList)
	# Take data set within specified offset
	offsetTimes <- times[startOffset:end]
	offsetData <- data[startOffset:end]
	inRangeTimes <- offsetTimes[1:end]
	inRangeData <- offsetData[1:end]

	# Mean of data
	meanPred <- myMean(inRangeData)

	# AR model order
	# AROrder <- 4
	for(i in seq(from=end, to=end, by=step)) {
		# Plot predicted data point for this time at previous fitting
		eval <- evalList[[i]]
		print(eval$optimRSquare)
		allEval <- eval$allEval$multiInf[1:end]
		print("allEval")
		ss <- ssError(allEval, inRangeData)
		rs <- rSquareError(allEval, inRangeData)
		rae <- rae(allEval, inRangeData)
		mad <- mad(allEval, inRangeData)
		mape <- mape(allEval, inRangeData)
		rmse <- rmse(allEval, inRangeData)
	}

	print(paste("SSE", ss))
	print(paste("RS", rs))	
	print(paste("RAE", rae))
	print(paste("MAD", mad))
	print(paste("MAPE", mape))
	print(paste("RMSE", rmse))

	# print(paste("EpiARSS", sseAR), quote=FALSE)
	# print(paste("ShiftSS", sseDiff), quote=FALSE)
	# print(paste("MeanARSS", sseMeanAR), quote=FALSE)

	# print(paste("EpiRS", rSqEpi), quote=FALSE)
	# print(paste("EpiARRS", rSqAR), quote=FALSE)
	# print(paste("ShiftRS", rSqDiff), quote=FALSE)
	# print(paste("MeanARRS", rSqMeanAR), quote=FALSE)
	
	# # print(paste("EpiARAIC", aicAR), quote=FALSE)
	# # print(paste("ShiftAIC", aicDiff), quote=FALSE)

	# print(paste("EpiARMAD", madAR), quote=FALSE)
	# print(paste("ShiftMAD", madDiff), quote=FALSE)

	# print(paste("EpiMean", myMean(abs(inRangeEvalPreds - inRangeData))))
	# print(paste("EpiARMean", myMean(abs(inRangeEvalPredsAR - inRangeData))))	
	# print(paste("ShiftMean", myMean(abs(inRangeShiftOffset - inRangeData))))

	# Set graph settings
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