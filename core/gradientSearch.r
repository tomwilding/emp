gradientSearch <- function(times, data) {
	# Fit linear gradient through series of datapoints
	intervalLength <- 5
	grad <- c()
	outbreakTimes <- c()
	k <- 10

	if (length(data) > intervalLength) {
		for (i in seq(from=1, to=(length(data) - intervalLength), by=1)) {
			interval <- times[i : (i + intervalLength)]
			intervalData <- data[i : (i + intervalLength)]
			linfit <- lm(intervalData ~ interval)
			# Store gradient change
			grad <- c(grad, unname(linfit$coef[2]))
		}
		# Find highest increases in gradients
		gradDiff <- grad[2:length(grad)] - grad[1:(length(grad) - 1)]
		orderedGradDiff <- sort(gradDiff, decreasing=TRUE)
		# Store outbreak times
		plot(times, data)
		for (i in 1:k) {
			maxDiffIndex <- which(gradDiff == orderedGradDiff[i])
			outbreakTimes <- c(outbreakTimes, maxDiffIndex)
			print(outbreakTimes)
		}
	}
}