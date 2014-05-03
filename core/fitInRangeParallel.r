fitInRangeParallel <- function(optimSIRMulti, i, offsetTimes, offsetData, initConds, initParams, epiTypes, k) {
	# Ensure length of range is greater than 0
	# assert_all_are_true(length(range) > 0)

	# Define eval vector
	eval <- c()

	# Truncate offsetData set
	truncTimes <- offsetTimes[1:i]
	truncData <- offsetData[1:i]

	# Fine Times for evaluation
	timeStep <- 0.05

	# Set optimParams to initParams to use if optimisation fails
	optimParams <- initParams
	# Range of feasible t0 values to explore referenced from offset offsetData
	# t0 <- proc.time()
	# For each possible start time optimise parameters of multiple epidemic - explore in parallel

	# EvaluateSIR over to plot fitted epidemics
	# lines(1:length(offsetTimes), offsetData, col='steelblue')
	# points(1:length(truncTimes), truncData, col='black', pch=16)
	###################################### Parallel evaluation at all feasible time points #######################################
	if (k > 1) {
		tryCatch({
			optimParams <- initParams
			# for (rep in 1 : 10) {
				optimParams <- optimSIRMulti(truncTimes, truncData, initConds, optimParams, epiTypes, k, i)
			# }
		}, warning = function(w) {
			print(w)
			print("optim warning")
		}, error = function(e) {
			print(e)
			print("optim failed")
		})
	}

	pastEval <- evalMulti(truncTimes, truncData, initConds, optimParams, epiTypes, k, 1, i)
	predInfectiousPast <- pastEval$multiInf
	# rSquare error to determine best time to start fitting
	rSquareError <- rSquareError(predInfectiousPast, truncData)

	# TODO: Don't want to restrict eval - eval over all points in evalMulti after optim over all but last n
	allEvalFine <- evalMulti(offsetTimes, offsetData, initConds, optimParams, epiTypes, k, timeStep, i)
	# Evaluate over all time
	allEval <- evalMulti(offsetTimes, offsetData, initConds, optimParams, epiTypes, k, 1, i) 
	startOffset <- offsets$startOffset
	endOffset <- offsets$endOffset

	# Plot inline for dev
 	fineTimes <- breakTime(offsetTimes, timeStep)
 	cl <- c("red","cyan","forestgreen","goldenrod2","red4")
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
 	for(k in 1:(length(allEvalFine$subInf))) {
 		sub <- allEvalFine$subInf[[k]]
 		subParams <- allEvalFine$subParams[[k]]
 		# Print sub epidemic graph
 		lines(fineTimes, sub, col=cl[k], lty=2)
 		lines(fineTimes, multiInf, col='black')
 		# lines(offsetTimes, multiInfCoarse, col='green')
 	}
 	dev.off()

	# Set values of eval
	eval$multiParams <- optimParams
	eval$initConds <- initConds
	# eval$optimTimes <- c(1, 23)
	eval$k <- k
	eval$optimRSquare <- rSquareError 
	# Get final residual from allEval infectious
	eval$finalResidual <- (offsetData[i] - allEval$multiInf[i])
	# Get all residuals up to current time
	eval$residuals <- (offsetData[1:i]) - (allEval$multiInf[1:i])
	# Set different multi and sub evals
	eval$nextPred <- allEval$multiInf[i+1]
	# eval$pastEval <- optimPastEval
	eval$allEval <- allEval
	eval$allEvalFine <- allEvalFine
	eval
}