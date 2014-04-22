findStartOffset <- function(data, minTruncation) {
	startOffset <- 1
	epidemicType <- 0
	epi <- FALSE
	lowerLimit <- 1000

	while(!epi) {
		epi <- TRUE
		meanData <- myMean(data[1:startOffset])
		sdData <- mySd(data[1:startOffset], meanData)
		lim <- 2 * sdData
		for (i in 1:minTruncation) {
			# Ensure all data points are above lim for epi
			epi <- epi && ((data[startOffset + i] - meanData) > lim)
		}
		startOffset <- startOffset + 1
	}
	startOffset - 1
}
