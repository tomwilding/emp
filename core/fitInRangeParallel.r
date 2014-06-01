fitInRangeParallel <- function(optimSIRMulti, i, offsetTimes, offsetData, initConds, initParams, epiTypes, ts, k, plotConfig, p) {

	# Define eval vector
	fit <- c()

	# Truncate offsetData set
	truncTimes <- offsetTimes[1:i]
	truncData <- offsetData[1:i]

	# Fine Times for evaluation
	timeStep <- 1
	fineTimeStep <- 0.05

	# Set optimParams to initParams to use if optimisation fails
	optimParams <- initParams

	if (k > 1) {
		# Find optimal parameters by optimising them to minimise the least square function
		tryCatch({
			optimParams <- optimSIRMulti(truncTimes, truncData, initConds, initParams, epiTypes, ts, k, fineTimeStep)
		}, warning = function(w) {
			print(w)
			print("optim warning")
		}, error = function(e) {
			print(e)
			print("optim failed")
		})
	}


	pastEval <- evalMulti(truncTimes, truncData, initConds, optimParams, epiTypes, ts, k, fineTimeStep)
	predInfectiousPast <- getObservations(pastEval$multiInf, fineTimeStep)
	# rSquare error to determine best time to start fitting
	rSquareError <- rSquareError(predInfectiousPast, truncData)
	# Evaluate over fine grain times
	allEvalFine <- evalMulti(offsetTimes, offsetData, initConds, optimParams, epiTypes, ts, k, fineTimeStep)
	# Get observations at observed data point times
	allEval <- getObservations(allEvalFine$multiInf, fineTimeStep)
	
	# Plot inline for dev
	if (p) {
	 	fineTimes <- breakTime(offsetTimes, fineTimeStep)
	 	cl <- c("red","cyan","forestgreen","goldenrod2","red4", "blue")
	 	setEPS()
	 	r <- plotConfig$run
	 	graphName <- paste("t", i, sep='')
	 	graphName <- paste(graphName, ".eps", sep='')
	 	postscript(paste(plotConfig$fileName, graphName, sep=''))	
	 	par(mar=c(7.1,4.1,4.1,2.1))
	 	plot(offsetTimes, offsetData, xlab='Epochs', ylab='Infected Individuals', col='steelblue')
	 	title(main=plotConfig$title, cex.main=1, cex.axis=0.8)
	 	daysText <- paste("Epochs after outbreak = ", i)
	 	mtext(daysText, 3, cex=0.8)
	 	# Plot offsetData points and actual offsetData lines
	 	lines(offsetTimes, offsetData, col='steelblue', lty=1)
	 	points(truncTimes, truncData, col='black', pch=16)
	 	# lines(fineTimes, allEvalFine$multiInf, lty=1)
	 	# multiInfCoarse <- allEval$multiInf
	 	multiInf <- allEvalFine$multiInf
	 	print("optimStartTimes")
	 	for(k in 1:(length(allEvalFine$subInf))) {
	 		sub <- allEvalFine$subInf[[k]]
	 		subParams <- allEvalFine$subParams[[k]]
	 		subStartTime <- allEvalFine$subStartTime[[k]]
	 		print(subStartTime)
	 		# Print sub epidemic graph
	 		lines(fineTimes, sub, col=cl[k], lty=2)
	 		lines(fineTimes, multiInf, col='black')
	 		# lines(offsetTimes, multiInfCoarse, col='green')
	 	}
	 	dev.off()
	}
	# Set values of fit
	fit$optimParams <- optimParams
	fit$k <- k
	fit$optimRSquare <- rSquareError 
	# Get final residual from allEval infectious
	fit$finalResidual <- (offsetData[i] - allEval[i])
	# Get all residuals up to current time
	fit$residuals <- (offsetData[1:i]) - (allEval[1:i])
	# Set different multi and sub fits
	fit$nextPred <- allEval[i+1]
	# fit$pastEval <- optimPastEval
	fit$allEval <- allEval
	fit$allEvalFine <- allEvalFine
	fit
}