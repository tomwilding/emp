# Helper function for LMS to compute sse
sseMulti <- function(params, times, data, initConds, epiTypes, ts, k) {
	# Calculate the sum of squared error at the current point
	outOfBounds <- FALSE
	subEpiNumParamsOffset <- 0
	
	for (i in 1:k) {
		# Get sub epidemic type and parameters
		subEpiNumParams <- epiTypes[i]
		paramsMulti <- params[(subEpiNumParamsOffset + 1) : (subEpiNumParamsOffset + subEpiNumParams)]
		initCondsMulti <- initConds[(subEpiNumParamsOffset + 1) : (subEpiNumParamsOffset + subEpiNumParams)]
		subEpiNumParamsOffset <- subEpiNumParamsOffset + subEpiNumParams
		if (subEpiNumParams == 1) {
			# Spike Epidemic
			# Set Spike epidemic optimised parameters
			gamma <- exp(paramsMulti[1])
			I0 <- exp(paramsMulti[2])
			if (gamma > 1 || gamma <= 1e-2) {
				sse <- Inf
				outOfBounds <- TRUE
			}
		}
		else if (subEpiNumParams == 3) {
			# Set SIR Epidemic optimised parameters
			beta <- exp(paramsMulti[1])
			gamma <- exp(paramsMulti[2])
			S0 <- exp(paramsMulti[3])
			I0 <- initCondsMulti[2]
			# Force optimisation to advance within parameter ranges
			if (beta > 1 || gamma > 1 || beta <= 1e-6 || gamma <= 1e-2 || S0 < I0) {
				sse <- Inf
				outOfBounds <- TRUE
			}
		}
	}
	if (!outOfBounds){
		granularity <- 1
		allPredInf <- evalMulti(times, data, initConds, params, epiTypes, ts, k, granularity)
		predInf <- allPredInf$multiInf
		sse <- sum((predInf - data[1:(length(data)-2)])^2)
	}
	sse
}