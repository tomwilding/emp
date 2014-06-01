evalMulti <- function(times, data, initConds, params, epiTypes, ts, k, timeStep) {
	require(deSolve)

	fineTimes <- breakTime(times, timeStep)
	predInfectious <- numeric(length(fineTimes))
	# Get initial I0
	I0 <- initConds[2]
	predI <- 0
	eval <- c()
	subEpiNumParamsOffset <- 0
	paramsMulti <- c()

	for (i in 1:k) {
		# Get sub epidemic type and parameters
		subEpiNumParams <- epiTypes[i]
		if (subEpiNumParams > 0) {
			paramsMulti <- params[(subEpiNumParamsOffset + 1) : (subEpiNumParamsOffset + subEpiNumParams)]
			initCondsMulti <- initConds[(subEpiNumParamsOffset + 1) : (subEpiNumParamsOffset + subEpiNumParams)]
			subEpiNumParamsOffset <- subEpiNumParamsOffset + subEpiNumParams

			# Evaluate epidemic according to type
			if (subEpiNumParams == 4) {
				epiStartTime <- logisticTransform(max(1,ts[i] - 20), paramsMulti[4], ts[i])
				# Update S0
				initCondsMulti[1] <- exp(paramsMulti[3])
				# Get predictions of SIR given current parameters
				preds <- ode(y=initCondsMulti, times=fineTimes, func=sir, parms=paramsMulti)
				predInf <- (preds[,3])
			} else if (subEpiNumParams == 1) {
				epiStartTime <- ts[i]
				# Update Spike epidemic parameters
				# Update I0 as unexplained prediction of Infectious for this epidemic
				initCondsMulti[1] <- I0
				preds <- ode(y=initCondsMulti, times=fineTimes, func=expDec, parms=paramsMulti)
				predInf <- preds[,2]
			}
		} else {
			# No epidemic detected
			epiStartTime <- 1
			predInf <- array(1, length(fineTimes)) * data[1]
		}

		# Find index to start fitting k+1 epidemic
		nearsetStartTime <- round(epiStartTime / timeStep) * timeStep
		t0Index <- ((nearsetStartTime - 1) / timeStep)
		zeros <- numeric(t0Index)
		predInf <- c(zeros, predInf)
		# Truncate to length of data
		predInf <- predInf[1:length(fineTimes)]
		# Add the predictions for the current epidemic to the previous predictions
		predInfectious <- predInfectious + predInf

		# Record sub epidemic parameters
		eval$subInf[[i]] <- predInf
		eval$subStartTime[[i]] <- epiStartTime
		eval$subParams[[i]] <- paramsMulti
		eval$multiInf <- predInfectious
		eval$multiParams <- params

		# Set I0 for next epidemic using combined predicted I0 at next t0
		if (i < k) {
			# Only predict I0 for kth epidemic
			t1Index <- which(fineTimes == times[ts[i + 1]])
			predI <- predInfectious[t1Index]
			# Update I0 as unexplained prediction of Infectious for the next epidemic
			I0 <- max(data[ts[i+1]] - predI, 1)
		}
	}
	eval
}