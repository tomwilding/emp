fitInRangeParallel <- function(optimSIRMulti, i, times, data, initConds, initParams, ts, k, range, actualFit, plotConfig) {
	require(doMC)
	require(assertive)
	cl <- c("red","cyan","forestgreen","goldenrod2","red4")
	# Ensure length of range is greater than 0
	assert_all_are_true(length(range) > 0)

	# Register multi core backend with n cores
	n <- length(range)
	registerDoMC(n)

	# Truncate data set
	truncTimes <- times[1:i]
	truncData <- data[1:i]

	# Fine Times for evaluation
	timeStep <- 0.05

	# Set optimParams to initParams to use if optimisation fails
	optimParams <- initParams
	# Range of feasible t0 values to explore referenced from offset data
	# t0 <- proc.time()
	# For each possible start time optimise parameters of multiple epidemic - explore in parallel

	# EvaluateSIR over to plot fitted epidemics
	lines(1:length(times), data, col='steelblue')
	points(1:length(truncTimes), truncData, col='black', pch=16)

	allTimeEval <- foreach (t=range) %dopar% {
		# time value t0 referenced from offset data
		tsExplore <- c(ts[1:k-1],t)
		# Find optimal beta and gamma by optimising them to minimise the least square function
		# OptimSIRMulti passed in from call to setSolver
		tryCatch({
			optimParams <- optimSIRMulti(truncTimes, truncData, initConds, initParams, tsExplore, k)
		}, warning = function(w) {
			print("optim warning")
		}, error = function(e) {
			print("optim failed")
		})
		allEval <- evalSIRMulti(truncTimes, truncData, initConds, optimParams, tsExplore, k, 1)
		predInfectiousPast <- allEval$multiInf 
		# rSquare error to determine best time to start fitting
		rSquareErrorPast <- rSquareError(predInfectiousPast, truncData)
		# Build list of optimisation results
		list(t, rSquareErrorPast, allEval)
	}

	# Get maximal rSquare index within parallel combined list
	maxRSIndex <- 1
	maxRS <- allTimeEval[[1]][[2]]
	for (i in 1:length(allTimeEval)) {
		RS <- allTimeEval[[i]][[2]]
		if (RS > maxRS) {
			maxRS <- RS
			maxRSIndex <- i
		}
	}

	optimTime <- allTimeEval[[maxRSIndex]][[1]]
	optimRSquare <- allTimeEval[[maxRSIndex]][[2]]
	# Optimal sub and combined epidemic parameters
	allEval <- allTimeEval[[maxRSIndex]][[3]]

	# Evaluate over fine granularity time
	optimParams <- allEval$multiParams
	allEvalFine <- evalSIRMulti(times, data, initConds, optimParams, c(ts[1:k-1], optimTime), k, timeStep)

	# setEPS()
	# graphName <- paste("t", i, sep='')
	# graphName <- paste(graphName, ".eps", sep='')
	# postscript(paste(plotConfig$fileName, graphName, sep=''))	
	# par(mar=c(7.1,4.1,4.1,2.1))
	# plot(1:length(times), data, xlab='Epochs', ylab='Infected Individuals', col='steelblue')
	# title(main=plotConfig$title, cex.main=1, cex.axis=0.8)
	# daysText <- paste("Epochs after outbreak = ", i)
	# mtext(daysText, 3, cex=0.8)
	# # lines(fineTimes, allEvalFine$multiInf, lty=1)
	# multiInf <- allEvalFine$multiInf
	# for(k in 1:(length(allEvalFine$subInf))) {
	# 	sub <- allEvalFine$subInf[[k]]
	# 	subParams <- allEvalFine$subParams[[k]]
	# 	# Print sub epidemic graph
	# 	lines(fineTimes, sub, col=cl[k], lty=2)
	# 	lines(fineTimes, multiInf, col="black")
	# }
	# dev.off()

	# allEval <- decomposeEpidemics(times, data, initConds, optimParams, c(ts[1:k-1], optimTime), k, actualFit, plotConfig)
	allEval$multiParams <- optimParams
	allEval$initConds <- initConds
	allEval$optimTime <- optimTime
	allEval$optimTimes <- c(ts[1:k-1], optimTime)
	allEval$k <- k
	allEval$optimRSquare <- optimRSquare
	allEval$allEvalFine <- allEvalFine
	allEval
}