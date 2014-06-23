reconstructPlot <- function(times, data, offsets, thresholds, initParams, initConds, plotConfig) {
	require("fanplot")
	require("forecast")
	# Unpack settings
	minTruncation <- offsets$minTruncation
	startOffset <- offsets$startOffset
	endOffset <- offsets$endOffset
	lim <- thresholds$lim
	
	# Colours for plotting
	cl <- c("red","deepskyblue3","forestgreen","goldenrod2","red4", "blue")

	# Loop through all objects

	end <- length(evalList)
	for(i in seq(from=minTruncation, to=end, by=1)) {
		# Set graph settings
		setEPS()
		graphName <- paste("t", i, sep='')
		graphName <- paste(graphName, ".eps", sep='')
		postscript(paste(plotConfig$fileName, graphName, sep=''))	
		
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
		epiTypes <- eval$epiTypes
		ts <- eval$optimTimes

		# Main plot
		par(mar=c(6.1,4.1,4.1,2.1))
		plot(offsetTimes, offsetData, xlab='Time (Days)', ylab='Infected Individuals', xaxs='i', col='steelblue')
		title(main=plotConfig$title, cex.main=0.9, cex.axis=0.8)
		# daysText <- paste("Epoch", i)
		# mtext(daysText, 3, cex=0.8)

		# Future RSqaure window to evaluate over
		NF <- 1
		
		rSquarePast <- eval$optimRSquare
		# rSquareNF <- rSquareError(allEval$multiInf[i:(i+NF)], offsetData[i:(i+NF)])
		if (!is.null(rSquarePast)) {
			RSqPastText <- paste("Past RSquare = ", signif(rSquarePast, digits=3))
			# RSqNFText <- paste("RSquare t+1 = ", signif(rSquareNF, digits=3))
			mtext(RSqPastText, 1, at=plotConfig$rat, padj=6, cex=0.7)
			# mtext(RSqNFText, 1, at=plotConfig$rat, padj=8, cex=0.7)
		}

		# Plot rectangle of past time first to plot on top
		rect(-1000000, -1000000, i + startOffset-1, 1000000, col=rgb(0.85,0.85,0.85), border=NA)
		# rect(i + startOffset-1, -1000000, min(i + startOffset-1 + NF, end + startOffset-1), 1000000, col=rgb(0.95,0.95,0.95), border=NA)
		# Plot all data points over rectangle
		points(offsetTimes, offsetData, col='steelblue')

		# Plot data points and actual data lines
		lines(offsetTimes, offsetData, col='steelblue', lty=1)
		
		# # Fan 
		# net <- ts(offsetData)
		# m <- auto.arima(net)
		# mm <- matrix(NA, nrow=1000, ncol=5)
		# 	for(mt in 1:1000)
  # 				mm[mt,] <- simulate(m, nsim=5)
		# fan(pn(mm), start=10, anchor=offsetData[10], type="interval", probs=seq(5, 95, 5), ln=c(50, 80))
		# abline(v=10)
		points(truncTimes, truncData, col='black', pch=16)
		# Plot lines using fine time for sub epidemics
		multiInf <- allEvalFine$multiInf
		for(k in 1:(length(allEvalFine$subInf))) {
			sub <- allEvalFine$subInf[[k]]
			subParams <- allEvalFine$subParams[[k]]
			# Print sub epidemic graph
			lines(fineTimes, sub, col=cl[k], lty=2)
			lines(fineTimes, multiInf, col="black")
			# Params - Add IO to S0
			if (k > 1) {
				epiType <- epiTypes[k]
				if (epiType == 3) {
					S0 <- exp(subParams[3])
					ParamText <- paste(c("Beta = ",", Gamma = ",", S0 = ", ", t0 = "), c(signif(exp(subParams[1:2]), digits=3), round(S0, digits=0), ts[k]), collapse='')
					mtext(ParamText, 1, at=plotConfig$pat, padj=4.5+(2*(k-2)), cex=0.7, col=cl[k])
				} else if (epiType == 1) {
					gamma <- exp(subParams[1])
					ParamText <- paste(c("Gamma = ", ", t0 = "), c(signif(exp(subParams[1]), digits=3), ts[k]), collapse='')
					mtext(ParamText, 1, at=plotConfig$pat, padj=4.5+(2*(k-2)), cex=0.7, col=cl[k])
				}
			}

			# Legend
			# legendText <- paste("Epidemic", 1:(k - 1))
			legendText <- c("Combined", "Baseline")
			lineType <- c(1,rep(2,k))
			col <- c(1,cl[1:k])
		}
		legend("topright",legendText, col=c("black", cl[1:(length(allEvalFine$subInf))]), lty=lineType, cex=0.8)
		dev.off()
	}
}