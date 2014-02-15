# Helper function for LMS to compute sse
sseSIRMulti <- function(params, times, data, initConds, ts, k) {
	# Calculate the sum of squared error at the current point
	inf <- 0
	for (i in 1:k) {
		paramsMulti <- params[(3*(i-1)+1):(3*i)]
		initCondsMulti <- initConds[(3*(i-1)+1):(3*i)]
		beta <- exp(paramsMulti[1])
		gamma <- exp(paramsMulti[2])
		S0 <- exp(paramsMulti[3])
		I0 <- initCondsMulti[2]
		if (beta>1 || gamma>1 || beta<=1e-6 || gamma<=1e-6 || S0<I0) {
			sse <- Inf
			inf <- 1
		}
	}
	if (inf==0){
		allPredInf <- evalSIRMulti(times, data, initConds, params, ts, k, 1)
		predInf <- allPredInf$multiInf
		sse <- sum((predInf - data)^2)
	}
	sse
}