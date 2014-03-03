evalSIRMulti <- function(times, data, initConds, params, epiTypes, ts, k, granularity) {
	require(deSolve)
	
	fineTimes <- breakTime(times, granularity);
	predInfectious <- numeric(length(fineTimes))
	I0 <- initConds[2]
	predI <- 0
	eval <- c()

	for (i in 1:k) {
		# Get sub epidemic type and parameters
		subEpiNumParams <- epiTypes[i]
		paramsMulti <- params[(subEpiNumParams * (i-1)+1) : (subEpiNumParams * i)]
		initCondsMulti <- initConds[(subEpiNumParams * (i-1)+1) : (subEpiNumParams * i)]
		if (subEpiNumParams > 2) {
			# Update SIR epidemic parameters
			# Update S0
			initCondsMulti[1] <- exp(paramsMulti[3])
			# Update I0
			initCondsMulti[2] <- I0
		}
		# Get predictions given current parameters
		preds <- as.data.frame(lsoda(y=initCondsMulti, times=fineTimes, func=sir, parms=paramsMulti))
		predInf <- (preds[,3])

		# If t0 > 1 then set offset in predicted infectious
		if (i > 1) {
			# Find index to start fitting k+1 epidemic
			t0Index <- which(fineTimes == (times[ts[i]]))
			# Offset
			zeros <- numeric(t0Index - 1)
			predInf <- c(zeros, predInf)
			# Truncate to length of data
			predInf <- predInf[1:length(fineTimes)]
		}
		# Add the predictions for the current epidemic to the previous predictions
		predInfectious <- predInfectious + predInf

		# Record sub epidemic parameters
		eval$subInf[[i]] <- predInf
		eval$subParams[[i]] <- c(paramsMulti, I0)
		eval$multiInf <- predInfectious
		eval$multiParams <- params
		# Set I0 for next epidemic using combined predicted I0 at next t0
		if (i < k) {
			# Check if S0 < I0 - why would I0 be less than S0
			# S0 optimised from optim so at the start S0 < I0 in.
			# Only predict I0 for kth epidemic
			t1Index <- which(fineTimes == times[ts[i+1]])
			predI <- predInfectious[t1Index]
			I0 <- max(data[ts[i+1]] - predI, 1)
		}
	}
	eval
}