gradientSearch <- function(times, data, plotConfig) {
	# Fit linear gradient through series of datapoints
	intervalLength <- 20
	grad <- c()
	outbreakTimes <- c()
	numSamples <- 5
	samplesTaken <- 0

 	setEPS()
 	graphName <- paste("gradSearch.eps")
 	postscript(paste(plotConfig$fileName, graphName, sep=''))	
	
 	truncData <- data[1:length(data)]
	truncTimes <- times[1:length(times)]
	
	if (length(truncData) > intervalLength) {
		for (i in seq(from=1, to=(length(truncData) - intervalLength), by=1)) {
			interval <- truncTimes[i : (i + intervalLength)]
			intervalData <- truncData[i : (i + intervalLength)]
			linfit <- lm(intervalData ~ interval)
			# Store gradient change
			grad <- c(grad, unname(linfit$coef[2]))
		}
		# Find highest increases from positive gradients
		gradDiff <- grad[2:length(grad)] - grad[1:(length(grad) - 1)]
		orderedGradDiff <- sort(gradDiff, decreasing=TRUE)
		# Store outbreak truncTimes
		plot(truncTimes, truncData, type='l')
		i <- 1
		while (samplesTaken < numSamples) {
			maxDiffIndex <- which(gradDiff == orderedGradDiff[i]) + intervalLength
			if (maxDiffIndex < length(grad)) {
				if (grad[maxDiffIndex] > 0 && newPoint(maxDiffIndex, outbreakTimes)) {
					outbreakTimes <- c(outbreakTimes, maxDiffIndex)
					samplesTaken <- samplesTaken + 1
					abline(v=maxDiffIndex)
				}
			}
			print(outbreakTimes)
			i <- i + 1
		}
	}
	title(main=plotConfig$title, cex.main=1, cex.axis=0.8)
	dev.off()
}

newPoint <- function(val, outbreakTimes) {
	minDistance <- 20
	new <- TRUE
	for (ot in outbreakTimes) {
		print(ot)
		new <- new && (abs(ot - val) > minDistance)
	}
	new
}
