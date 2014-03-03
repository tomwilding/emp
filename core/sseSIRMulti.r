# Helper function for LMS to compute sse
sseSIRMulti <- function(params, times, data, initConds, epiTypes, ts, k) {
	# Calculate the sum of squared error at the current point
	inf <- 0
	for (i in 1:k) {
		# Get sub epidemic type and parameters
		subEpiNumParams <- epiTypes[i]
		paramsMulti <- params[(subEpiNumParams * (i-1)+1) : (subEpiNumParams * i)]
		initCondsMulti <- initConds[(subEpiNumParams * (i-1)+1) : (subEpiNumParams * i)]
		if (subEpiNumParams == 2) {
			# Spike Epidemic
			# TODO: Set Spike epidemic optimised parameters
		}
		else if (subEpiNumParams > 2) {
			# Set SIR Epidemic optimised parameters
			beta <- exp(paramsMulti[1])
			gamma <- exp(paramsMulti[2])
			S0 <- exp(paramsMulti[3])
			I0 <- initCondsMulti[2]
			# Force optimisation to advance within parameter ranges
			if (beta>1 || gamma>1 || beta<=1e-6 || gamma<=1e-6 || S0<I0) {
				sse <- Inf
				inf <- 1
			}
		}
	}
	if (inf==0){
		granularity <- 1
		allPredInf <- evalSIRMulti(times, data, initConds, params, epiTypes, ts, k, granularity)
		predInf <- allPredInf$multiInf
		sse <- sum((predInf - data)^2)
	}
	sse
}