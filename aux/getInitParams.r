getInitParams <- function(optimMethod, i, offsetTimes, offsetData, initConds, initParams, epiTypes, ts, k, plotConfig, data) {
	# Fit over range of new params and take best
	numIntervals <- 10
    minT <- 0; maxT <- 6
    tValsPos <- seq(from=minT, to=maxT, by=(maxT - minT) / numIntervals)
    tVals <- tValsPos - ((maxT - minT) / 2)
	maxRS <- 0
	# Remove last epidemic parameters
	initParams <- reduceParams(initParams, epiTypes[length(epiTypes)])
	print("Finding new initParams")
	for (t in tVals) {
		print(t)
		newParams <- c(initParams, c(log(0.001), log(0.01), log(1000), t))
		print(newParams)
		eval <- fitInRangeParallel(setSolver("FS", k, epiTypes), i, offsetTimes, offsetData, initConds, newParams, epiTypes, ts, k, plotConfig, 1)
		print("afterEval")
		if (eval$optimRSquare > maxRS) {
			maxRS <- eval$optimRSquare
			initParamsMax <- newParams
		}
		newParams <- c()
	}
	print("maxInitParams")
	print(initParamsMax)
	readline()
	# readline()
	initParamsMax
}