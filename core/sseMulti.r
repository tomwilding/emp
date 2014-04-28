# Helper function for LMS to compute sse
sseMulti <- function(params, times, data, initConds, epiTypes, k, tmax) {
	# Calculate the sum of squared error at the current point
	outOfBounds <- FALSE
	subEpiNumParamsOffset <- 0
	
	for (i in 1:k) {
		# Get sub epidemic type and parameters
		subEpiNumParams <- epiTypes[i]
		paramsMulti <- params[(subEpiNumParamsOffset + 1) : (subEpiNumParamsOffset + subEpiNumParams)]
		initCondsMulti <- initConds[(subEpiNumParamsOffset + 1) : (subEpiNumParamsOffset + subEpiNumParams)]
		subEpiNumParamsOffset <- subEpiNumParamsOffset + subEpiNumParams

		# Check parameters are in bounds
		if (subEpiNumParams == 1) {
			# Spike Epidemic
			# Set Spike epidemic optimised parameters
			gamma <- exp(paramsMulti[1])
			I0 <- exp(paramsMulti[2])
			if (gamma > 1 || gamma <= 1e-6) {
				sse <- Inf
				outOfBounds <- TRUE
			}
		}
		else if (subEpiNumParams == 4) {
			# Set SIR Epidemic optimised parameters
			beta <- exp(paramsMulti[1])
			gamma <- exp(paramsMulti[2])
			S0 <- exp(paramsMulti[3])
			# t0 <- logisticTransform(paramsMulti[4], tmax)
			I0 <- initCondsMulti[2]
			# Force optimisation to advance within parameter ranges
			if (beta > 1 || gamma > 1 || beta <= 1e-6 || gamma <= 1e-6 || S0 < I0) {
				sse <- Inf
				outOfBounds <- TRUE
			}
		}
	}
	# If parameters are in bounds then eval and get sse
	if (!outOfBounds){
		granularity <- 1
		eval <- evalMulti(times, data, initConds, params, epiTypes, k, granularity, tmax)
		predInf <- eval$multiInf
		sse <- sum((predInf - data)^2)
	}
	sse
}