fitInRangeParallel <- function(optimSIRMulti, i, times, data, initConds, initParams, ts, k, range, actualFit, plotConfig) {
	require(doMC)
	require(assertive)

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
	fineTimes <- breakTime(times, timeStep);

	# Set optimParams to initParams to use if optimisation fails
	optimParams <- initParams
	# Range of feasible t0 values to explore referenced from offset data
	# t0 <- proc.time()
	allOptimParams <- foreach (t=range, .combine="rbind") %dopar% {
		# time value t0 referenced from offset data
		tsExplore <- c(ts[1:k-1],t)
		# Find optimal beta and gamma by optimising them to minimise the least square function - OptimSIRMulti passed in from call to setSolver
		tryCatch({
			optimParams <- optimSIRMulti(truncTimes, truncData, initConds, initParams, tsExplore, k)
		}, warning = function(w) {
			print("optim warning")
		}, error = function(e) {
			print("optim failed")
		})
		predInfectiousPast <- evalSIRMulti(truncTimes, truncData, initConds, optimParams, tsExplore, k)
		# rSquare error to determine best time to start fitting
		rSquareErrorPast <- rSquareError(predInfectiousPast, truncData)
		# Build list of optimisation results
		c(t, rSquareErrorPast, optimParams)
	}

	# Find max rSqaureError entry
	if (length(range) > 1) {
		# Find index with best fit according to max RSq
		optimIndex <- which.max(allOptimParams[,2])
		# Get max parameters
		optimList <- allOptimParams[optimIndex,]
	} else {
		# Base case when vector of single results
		optimList <- allOptimParams
	}

	optimTime <- optimList[1]
	optimRSquare <- optimList[2]
	# Index optim params to length of list to get all params
	optimParams <- optimList[3:length(optimList)]

	# EvaluateSIR over to plot fitted epidemics
	par(mar=c(7.1,4.1,4.1,2.1))
	plot(1:length(times), data, xlab='Epochs', ylab='Infected Individuals', col='steelblue')
	title(main=plotConfig$title, cex.main=1, cex.axis=0.8)
	daysText <- paste("Epochs after outbreak = ", i)
	mtext(daysText, 3, cex=0.8)
	RSqPastText <- paste("Past RSquare = ", signif(optimRSquare, digits=3))
	mtext(RSqPastText, 1, at=plotConfig$rat, padj=6, cex=0.7) #other
	lines(1:length(times), data, col='steelblue')
	points(1:length(truncTimes), truncData, col='black', pch=16)

	allPredInf <- decomposeEpidemics(times, data, initConds, optimParams, c(ts[1:k-1],optimTime), k, actualFit, plotConfig)
	allPredInf$optimTime <- optimTime
	allPredInf$optimParams <- optimParams
	allPredInf$optimRSquare <- optimRSquare
	allPredInf
}