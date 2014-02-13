decomposeEpidemics <- function(times, data, initConds, params, ts, k, actualFit, plotConfig) {
	# source('sir.r')
	# Colours for plotting
	cl <- c("red","cyan","forestgreen","goldenrod2","red4")
	# Fine Times for evaluation
	timeStep <- 0.05
	fineTimes <- breakTime(times, timeStep);
	fineTimesPlot <- breakTime(c(1:length(times)), timeStep)
	predInfectious <- numeric(length(fineTimes))
	I0 <- initConds[2]
	predI <- 0
	allPredInf <- c()

	for (i in 1:k) {
		paramsMulti <- params[(3*(i-1)+1):(3*i)]
		initCondsMulti <- initConds[(3*(i-1)+1):(3*i)]
		# Update S0
		initCondsMulti[1] <- exp(paramsMulti[3])
		# Update I0
		initCondsMulti[2] <- I0
		# Get predictions given current parameters
		preds <- as.data.frame(lsoda(y=initCondsMulti, times=fineTimes, sir, parms=paramsMulti))
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

		plotParams <- c(signif(exp(paramsMulti[1:2]), digits=3), round(exp(paramsMulti[3])+I0, , digits=0))
		ParamText <- paste(c("Beta = ",", Gamma = ",", S0 = "), plotParams, collapse='')
		mtext(ParamText, 1, at=plotConfig$pat, padj=6+(2*(i-1)), cex=0.7, col=cl[i])

		# Print update
		if (actualFit == 1) {
			lines(fineTimesPlot, predInf, col=cl[i], lty=2)
			allPredInf$sub[[i]] <- predInf
			allPredInf$plotParams[[i]] <- plotParams
			# allPredInf$optimParams <- c(paramsMulti[1:2], paramsMulti[3]+I0)	
		}

		# Set I0 for next epidemic using combined predicted I0 at next t0
		if (i < k) {
			t1Index <- which(fineTimes == times[ts[i+1]])
			predI <- predInfectious[t1Index]
			I0 <- max(data[ts[i+1]] - predI, 1)
		}
	}
	if (actualFit == 1) {
		lines(fineTimesPlot, predInfectious, lty=actualFit)
		allPredInf$overall <- predInfectious
	}
	legendText <- paste("Epidemic",1:k)
	legendText <- c("Combined", legendText)
	lineType <- c(1,rep(2,k))
	legend("topright",legendText, col=c(1,cl[1:k]), lty=lineType, cex=0.8)
	dev.off()
	allPredInf
}