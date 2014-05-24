getInitParams <- function(optimMethod, i, offsetTimes, offsetData, initConds, initParams, epiTypes, ts, k, plotConfig, data) {
	# Fit over range of new params and take best
	minB <- 1e-5; maxB <- 1
	minG <- 1e-2; maxG <- 1
	betaVals <- seq(from=minB + ((maxB - minB) / 4), to=maxB - ((maxB - minB) / 4), by=((maxB - minB) / 4))
	gammaVals <- seq(from=minG + ((maxG - minG) / 4), to=maxG - ((maxG - minG) / 4), by=((maxG - minG) / 4))
	initParamsMax <- c(initParams, c(minB, minG, incOrderOf(data[i]), 0))
	maxRS <- 0
	print("Finding new initParams")
	for (b in betaVals) {
		for (g in gammaVals) {
			newParams <- c(initParams, c(logit(minB, b, maxB), logit(minG, g, maxG), logit(1e0, incOrderOf(data[i]), 1e6), 0))
			eval <- fitInRangeParallel(optimMethod, i, offsetTimes, offsetData, initConds, newParams, epiTypes, ts, k, plotConfig, 0)
			if (eval$optimRSquare > maxRS) {
				maxRS <- eval$optimRSquare
				initParamsMax <- newParams
			}
		}
	}
	print("maxInitParams")
	print(initParamsMax)
	initParamsMax
}