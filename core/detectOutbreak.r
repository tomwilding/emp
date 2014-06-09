detectOutbreak <- function(residuals, nRes, startTime, k) {
	outbreak <- 0
	# Ensure more than one residual before the last n residuals to calculate sdRes
	resLength <- length(residuals)
	if (resLength > nRes + 1) {
		# Get standard deviation of residuals before the ones considered
		inRangeResiduals <- residuals[startTime : (resLength - nRes)]
		meanRes <- myMean(inRangeResiduals)
		sdRes <- mySd(inRangeResiduals, meanRes)

		print(paste("meanRes", meanRes))
		print(paste("sdRes ", sdRes))
		print(paste("Outbreaklim", meanRes + sdRes * 2))
		print(paste("Explim", meanRes + sdRes * 6))
		# If current residual sd is above zero check if last n residuals are above set number of sd
		if (sdRes > 0) {
			# Index of first residual to check
			startResIndex <- resLength - nRes + 1
			# Assume incRes is True and check condition for all n residuals
			outbreakRes <- min(residuals[startResIndex : resLength])
			expRes <- residuals[resLength]
			print(paste("OutbreakRes", outbreakRes))
			print(paste("ExpRes", expRes))
			# Set epidemic type according to residual limit
			outbreakLim <- (meanRes + (sdRes * 2))
			expLim <- (meanRes + (sdRes * 6))
			# If minimum residual increase is more than required, then set type
			if (outbreakRes > outbreakLim) {
				outbreak <- 3
			} else if (expRes > expLim) {
				outbreak <- 1
			}
		}
	}
	outbreak	
}