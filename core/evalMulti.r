evalMulti <- function(times, data, initConds, params, epiTypes, ts, k, granularity) {
	require(deSolve)

	fineTimes <- breakTime(times, granularity);
	predInfectious <- numeric(length(fineTimes))
	# Get initial I0
	# Update for next epidemic according to unexplained offset from previous epidemic
	I0 <- initConds[2]
	predI <- 0
	eval <- c()
	eval$subParams[[1]] <- c(1,1,1)
	subEpiNumParamsOffset <- 0
	paramsMulti <- c()

	for (i in 1:k) {
		# Get sub epidemic type and parameters
		subEpiNumParams <- epiTypes[i]
		paramsMulti <- params[(subEpiNumParamsOffset + 1) : (subEpiNumParamsOffset + subEpiNumParams)]
		initCondsMulti <- initConds[(subEpiNumParamsOffset + 1) : (subEpiNumParamsOffset + subEpiNumParams)]
		subEpiNumParamsOffset <- subEpiNumParamsOffset + subEpiNumParams

		if (i > 1) {
			# Evaluate epidemic according to type
			if (subEpiNumParams == 3) {
				# Update SIR epidemic parameters
				# Update S0
				initCondsMulti[1] <- exp(paramsMulti[3])
				# print(initCondsMulti)
				# Update I0 computed using previous sub epidemics
				# initCondsMulti[2] <- I0
				# Get predictions of SIR given current parameters
				preds <- lsoda(y=initCondsMulti, times=fineTimes, func=sir, parms=paramsMulti)
				predInf <- (preds[,3])
			} else if (subEpiNumParams == 1) {
				# Update Spike epidemic parameters
				# Update I0 as unexplained prediction of Infectious for this epidemic
				initCondsMulti[1] <- I0
				preds <- lsoda(y=initCondsMulti, times=fineTimes, func=expDec, parms=paramsMulti)
				predInf <- preds[,2]
				# Update I0 in conds
				initConds[subEpiNumParamsOffset + 1] <- I0
			}
		} else {
			# Optimise baseline
			baseline <- exp(paramsMulti[1])
			predInf <- array(baseline, length(fineTimes))
		}
		optimConds <- initConds

		# Find index to start fitting k+1 epidemic
		t0Index <- which(fineTimes == (times[ts[i]]))
		# Offset
		zeros <- numeric(t0Index - 1)
		predInf <- c(zeros, predInf)
		# Truncate to length of data
		predInf <- predInf[1:length(fineTimes)]
		# Add the predictions for the current epidemic to the previous predictions
		predInfectious <- predInfectious + predInf

		# Record sub epidemic parameters
		eval$subInf[[i]] <- predInf
		eval$subParams[[i]] <- paramsMulti
		eval$multiInf <- predInfectious
		eval$multiParams <- params
		eval$multiConds <- optimConds
		# Set I0 for next epidemic using combined predicted I0 at next t0
		# TODO: Update I0 after or before setting it in eval??
		if (i < k) {
			# TODO: I0 only calculated from previous epidemic
			# Only predict I0 for kth epidemic
			t1Index <- which(fineTimes == times[ts[i + 1]])
			predI <- predInfectious[t1Index]
			# Update I0 as unexplained prediction of Infectious for the next epidemic
			I0 <- max(data[ts[i+1]] - predI, 1)
		}
	}
	eval
}