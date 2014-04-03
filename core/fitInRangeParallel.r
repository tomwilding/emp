fitInRangeParallel <- function(optimSIRMulti, i, offsetTimes, offsetData, initConds, initParams, epiTypes, ts, k, range, plotConfig) {
	require(doMC)

	# Ensure length of range is greater than 0
	# assert_all_are_true(length(range) > 0)

	# Define eval vector
	eval <- c()

	# Register multi core backend with n cores
	n <- length(range)
	registerDoMC(n)

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
	EvalOverTime <- foreach (t=range) %dopar% {
		# time value t0 referenced from offset offsetData
		tsExplore <- c(ts[1:k-1],t)
		# Find optimal beta and gamma by optimising them to minimise the least square function
		# OptimSIRMulti passed in from call to setSolver
		tryCatch({
			optimParams <- optimSIRMulti(truncTimes, truncData, initConds, initParams, epiTypes, tsExplore, k)
		}, warning = function(w) {
			print(w)
			print("optim warning")
		}, error = function(e) {
			print(e)
			print("optim failed")
		})
		pastEval <- evalMulti(truncTimes, truncData, initConds, optimParams, epiTypes, tsExplore, k, 1)
		predInfectiousPast <- pastEval$multiInf
		# rSquare error to determine best time to start fitting
		rSquareErrorPast <- rSquareError(predInfectiousPast, truncData)
		# Build list of optimisation results
		list(t, rSquareErrorPast, pastEval)
	}

	# Get maximal rSquare index within parallel combined list
	maxRSIndex <- 1
	maxRS <- EvalOverTime[[1]][[2]]
	for (r in 1:length(EvalOverTime)) {
		RS <- EvalOverTime[[r]][[2]]
		if (RS > maxRS) {
			maxRS <- RS
			maxRSIndex <- r
		}
	}

	# Get optimal values stored in parallel evaluation loop
	optimTime <- EvalOverTime[[maxRSIndex]][[1]]
	optimRSquare <- EvalOverTime[[maxRSIndex]][[2]]
	# Optimal sub and combined epidemic parameters
	optimPastEval <- EvalOverTime[[maxRSIndex]][[3]]
	# Evaluate over all fine granularity time
	optimParams <- optimPastEval$multiParams
	# TODO: Don't want to restrict eval - eval over all points in evalMulti after optim over all but last n
	allEvalFine <- evalMulti(offsetTimes, offsetData, initConds, optimParams, epiTypes, c(ts[1:k-1], optimTime), k, timeStep)
	# Evaluate over all time
	allEval <- evalMulti(offsetTimes, offsetData, initConds, optimParams, epiTypes, c(ts[1:k-1], optimTime), k, 1) 
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
 	plot(1:length(offsetTimes), offsetData, xlab='Epochs', ylab='Infected Individuals', col='steelblue')
 	title(main=plotConfig$title, cex.main=1, cex.axis=0.8)
 	daysText <- paste("Epochs after outbreak = ", i)
 	mtext(daysText, 3, cex=0.8)
 	# Plot offsetData points and actual offsetData lines
 	lines(offsetTimes, offsetData, col='steelblue', lty=1)
 	points(truncTimes, truncData, col='black', pch=16)
 	# lines(fineTimes, allEvalFine$multiInf, lty=1)
 	multiInf <- allEvalFine$multiInf
 	for(k in 1:(length(allEvalFine$subInf))) {
 		sub <- allEvalFine$subInf[[k]]
 		subParams <- allEvalFine$subParams[[k]]
 		# Print sub epidemic graph
 		lines(fineTimes, sub, col=cl[k], lty=2)
 		lines(fineTimes, multiInf, col="black")
 	}
 	dev.off()

	# Set values of eval
	eval$multiParams <- optimParams
	eval$initConds <- initConds
	eval$optimTime <- optimTime
	eval$optimTimes <- c(ts[1:k-1], optimTime)
	eval$k <- k
	eval$optimRSquare <- optimRSquare 
	# Get final residual from allEval infectious
	eval$finalResidual <- (offsetData[i] - allEval$multiInf[i])
	# Get all residuals up to current time
	eval$residuals <- (offsetData[1:i]) - (allEval$multiInf[1:i])
	# Set different multi and sub evals
	eval$pastEval <- optimPastEval
	eval$allEval <- allEval
	eval$allEvalFine <- allEvalFine
	eval
}