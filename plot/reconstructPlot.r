reconstructPlot <- function(times, data, offsets, thresholds, initParams, initConds, plotConfig) {
	# Unpack settings
	minTruncation <- offsets$minTruncation
	startOffset <- offsets$startOffset
	endOffset <- offsets$endOffset
	lim <- thresholds$lim
	
	# Colours for plotting
	cl <- c("red","cyan","forestgreen","goldenrod2","red4", "blue")

	# Loop through all objects

	end <- length(evalList)
	for(i in minTruncation:end) {
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
		# Print main graph

		# Main plot
		par(mar=c(6.1,4.1,4.1,2.1))
		plot(offsetTimes, offsetData, xlab='Epochs', ylab='Infected Individuals', col='steelblue')
		title(main=plotConfig$title, cex.main=0.9, cex.axis=0.8)
		daysText <- paste("Epochs after outbreak = ", i)
		mtext(daysText, 3, cex=0.8)

		# RSquare and labels
		# Future RSqaure window to evaluate over
		NF <- 7
		# Calculate RSqaure
		rSquarePast <- eval$optimRSquare
		
		# rSquareNF <- rSquareError(allEval$multiInf[i:(i+NF)], offsetData[i:(i+NF)])
		# RSqPastText <- paste("RSquare Past = ", signif(rSquarePast, digits=3))
		# RSqT1Text <- paste("RSquare Near Future = ", signif(rSquareNF, digits=3))
		# mtext(RSqPastText, 1, at=plotConfig$rat, padj=6, cex=0.7)
		# mtext(RSqT1Text, 1, at=plotConfig$rat, padj=8, cex=0.7)

		# Plot rectangle of past time first to plot on top
		rect(-1000000, -1000000, i + startOffset-1, 1000000, col=rgb(0.8,0.8,0.8), border=NA)
		rect(i + startOffset-1, -1000000, min(i + startOffset-1 + NF, end + startOffset-1), 1000000, col=rgb(0.95,0.95,0.95), border=NA)
		# Plot all data points over rectangle
		points(offsetTimes, offsetData, col='steelblue')

		# Plot data points and actual data lines
		lines(offsetTimes, offsetData, col='steelblue', lty=1)
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
			# S0 <- exp(subParams[3]) + subParams[4]
			# ParamText <- paste(c("Beta = ",", Gamma = ",", S0 = "), c(signif(exp(subParams[1:2]), digits=3), round(S0, digits=0)), collapse='')
			# mtext(ParamText, 1, at=plotConfig$pat, padj=6+(2*(k-1)), cex=0.7, col=cl[k])		

			# Legend
			legendText <- paste("Epidemic",1:k)
			legendText <- c("Combined", legendText)
			lineType <- c(1,rep(2,k))
			col <- c(1,cl[1:k])
		}
		legend("topright",legendText, col=c(1,cl[1:(length(allEvalFine$sub))]), lty=lineType, cex=0.8)
		dev.off()
	}
}