getInitParams <- function(optimMethod, i, offsetTimes, offsetData, initConds, initParams, epiTypes, ts, k, plotConfig, data) {
	# Fit over range of new params and take best
	minB <- 1e-5; maxB <- 1
	minG <- 1e-2; maxG <- 1
	minS <- 1e1; maxS <- 1e5
	numIncrements <- 4
	SVals <- seq(from=minS, to=maxS, by=((maxS - minS) / numIncrements))
	betaVals <- seq(from=minB, to=maxB, by=((maxB - minB) / numIncrements))
	gammaVals <- seq(from=minG, to=maxG, by=((maxG - minG) / numIncrements))
	initParamsMax <- c(initParams, c(minB, minG, minS, 0))
	maxRS <- 0
	print("Finding new initParams")
	for (s in SVals) {
		print(s)
		for (b in betaVals) {
			for (g in gammaVals) {
				newParams <- c(initParams, c(log(b), log(g), log(s), logit(max((i - 10), 1), (i - minTruncation), i)))
				eval <- fitInRangeParallel(optimMethod, i, offsetTimes, offsetData, initConds, newParams, epiTypes, ts, k, plotConfig, 0)
				if (eval$optimRSquare > maxRS) {
					maxRS <- eval$optimRSquare
					initParamsMax <- newParams
				}
			}
		}
	}
	print("maxInitParams")
	print(initParamsMax)
	# readline()
	initParamsMax
}