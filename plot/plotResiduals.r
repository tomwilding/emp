plotResiduals <- function(times, data, offset, thresholds, initParams, initConds, plotConfig) {
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
		
		# Graph settings
		setEPS()
		
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
		res <- eval$residuals
		if (i==length(evalList)) {
			postscript(paste(plotConfig$fileName, "acf"))
			acf(res)
			dev.off()

			postscript(paste(plotConfig$fileName, "pacf"))
			pacf(res)
			dev.off()

			# postscript(paste(plotConfig$fileName, "ar"))
			# ars <- ar(res)
			# print(ars)
			# future <- predict(ars, n.ahead=250)
			# all <- c(res, as.vector(future$pred))
			# plot(c(1:length(all)), all, type="l")
			# title("Autoregression")
			# dev.off()

			# Save res to csv
			write.table(file="finalRes.csv", res, quote=FALSE, sep=",")

			# Plot res
			postscript(paste(plotConfig$fileName, "res"))
			plot(1:length(res), res, type="l")
			dev.off()
		}

		# graphName <- paste("res", i, sep='')
		# graphName <- paste(graphName, ".eps", sep='')
		# postscript(paste(plotConfig$fileName, graphName, sep=''))	
		
		# # Plot the graph
		# plot(truncTimes, res, type="l") 

		# dev.off()
	}

}