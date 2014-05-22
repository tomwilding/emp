breakTime <- function(times, timeStep) {
	minTime <- min(times)
	maxTime <- max(times)
	fineTimes <- seq(minTime, maxTime, by=timeStep)
}