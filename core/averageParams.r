averageParams <- function(initParams, optimParamsList, epiTypesList, i, k) {
	# Find a new set of parameter seed using previous optimal parameter values
	# Average all epidmeic parameters independently
	averageB <- 0
	averageG <- 0
	averageS <- 0

	totalB <- 0
	totalG <- 0
	totalS <- 0

	subEpiNumParamsOffset <- 0
	window <- 10

	if (k > 1) {
		if (i > window) {
			for (p in (i - window - 1) : i) {
				optimParams <- optimParamsList[[i - 1]]
				epiTypes <- epiTypesList[[i - 1]]
				for (t in epiTypes) {
					if (t == 3) {
						totalB <- totalB + optimParams[subEpiNumParamsOffset + 1]
						totalG <- totalG + optimParams[subEpiNumParamsOffset + 2]
						totalS <- totalS + optimParams[subEpiNumParamsOffset + 3]
					}
					subEpiNumParamsOffset <- subEpiNumParamsOffset + t
				}
			}
			averageB <- totalB / window
			averageG <- totalG / window
			averageS <- totalS / window
		}
		# Need to average each epidemics parameters individually
		# averageParams <- 
	} else {
		averageParams <- initParams
	}
	averageParams
}