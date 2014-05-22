getObservations <- function(data, granularity) {
	observationPoints <- c()
	for (i in 1:((length(data)-1) / (1/granularity) + 1)) {
		observationPoints[i] <- data[(i - 1) * (1/granularity) + 1]
	}
	observationPoints
}