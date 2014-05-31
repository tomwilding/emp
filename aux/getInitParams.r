getInitParams <- function(optimMethod, i, offsetTimes, offsetData, initConds, initParams, epiTypes, ts, k, plotConfig, data) {
	# Fit over range of new params and take best
	minB <- 1e-6; maxB <- 1
	minG <- 1e-5; maxG <- 1
	minS <- 1e1; maxS <- 1e4
	minT <- 0; maxT <- 200;
	numIncrements <- 4
	SVals <- seq(from=minS, to=maxS, by=((maxS - minS) / numIncrements))
	betaVals <- seq(from=minB, to=maxB, by=((maxB - minB) / numIncrements))
	gammaVals <- seq(from=minG, to=maxG, by=((maxG - minG) / numIncrements))
	tVals <- seq(from=minT, to=maxT, by=((maxT - minT) / numIncrements)) - 100
	newParams <- c()
	maxRS <- 0
	print("Finding new initParams")
	inc <- 0
	for (t in tVals) {
		print(inc)
		for (s in SVals) {
			for (g in gammaVals) {
				for (b in betaVals) {
					for (r in 2:k) {
						newParams <- c(newParams, c(log(b), log(g), log(s), t))
					}
					eval <- fitInRangeParallel(setSolver("SP", k, epiTypes), i, offsetTimes, offsetData, initConds, newParams, epiTypes, ts, k, plotConfig, 0)
					print("eval")
					if (eval$optimRSquare > maxRS) {
						maxRS <- eval$optimRSquare
						initParamsMax <- newParams
					}
					newParams <- c()
				}
			}
		}
		inc <- inc + 1
	}
	print("maxInitParams")
	print(initParamsMax)
	# readline()
	initParamsMax
}