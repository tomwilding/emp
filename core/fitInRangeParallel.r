fitInRangeParallel <- function(optimSIRMulti, i, offsetTimes, offsetData, initConds, initParams, epiTypes, ts, k, range, plotConfig, p) {
	require(doMC)

	# Ensure length of range is greater than 0
	# assert_all_are_true(length(range) > 0)

	# Define eval vector
	eval <- c()

	# Register multi core backend with n cores
	registerDoMC(min(length(range), 60))

	# Truncate offsetData set
	truncTimes <- offsetTimes[1:i]
	truncData <- offsetData[1:i]

	# Fine Times for evaluation
	timeStep <- 0.05

	# Set optimParams to initParams to use if optimisation fails
	optimParams <- initParams

	###################################### Parallel evaluation at all feasible time points #######################################
	evalOverTime <- foreach (t=range) %dopar% {
		tsExplore <- c(ts[1:(k - 1)],t)
		# time value t0 referenced from offset offsetData
		# if (k > 1) {
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
		# }
		pastEval <- evalMulti(truncTimes, truncData, initConds, optimParams, epiTypes, tsExplore, k, 1)
		predInfectiousPast <- pastEval$multiInf
		# rSquare error to determine best time to start fitting
		rSquareErrorPast <- rSquareError(predInfectiousPast, truncData)
		# Build list of optimisation results
		list(t, rSquareErrorPast, pastEval)
	}

	# Get maximal rSquare index within parallel combined list
	maxRSIndex <- 1
	maxRS <- evalOverTime[[1]][[2]]
	for (r in 1:length(evalOverTime)) {
		RS <- evalOverTime[[r]][[2]]
		if (RS > maxRS) {
			maxRS <- RS
			maxRSIndex <- r
		}
	}

	# Get optimal values stored in parallel evaluation loop
	optimTime <- evalOverTime[[maxRSIndex]][[1]]
	optimRSquare <- evalOverTime[[maxRSIndex]][[2]]
	# Optimal sub and combined epidemic parameters
	optimPastEval <- evalOverTime[[maxRSIndex]][[3]]
	# Evaluate over all fine granularity time
	optimParams <- optimPastEval$multiParams
	optimConds <- optimPastEval$multiConds
	# TODO: Don't want to restrict eval - eval over all points in evalMulti after optim over all but last n
	allEvalFine <- evalMulti(offsetTimes, offsetData, initConds, optimParams, epiTypes, c(ts[1:(k-1)], optimTime), k, timeStep)
	# Evaluate over all time
	allEval <- evalMulti(offsetTimes, offsetData, initConds, optimParams, epiTypes, c(ts[1:(k-1)], optimTime), k, 1) 

	# Plot inline for dev
	# if (p) {
	#  	fineTimes <- breakTime(offsetTimes, timeStep)
	#  	cl <- c("red","cyan","forestgreen","goldenrod2","red4")
	#  	setEPS()
	#  	r <- plotConfig$run
	#  	graphName <- paste("t", i, sep='')
	#  	graphName <- paste(graphName, ".eps", sep='')
	#  	postscript(paste(plotConfig$fileName, graphName, sep=''))	
	#  	par(mar=c(7.1,4.1,4.1,2.1))
	#  	plot(offsetTimes, offsetData, xlab='Days', ylab='Infected Individuals', col='steelblue')
	#  	title(main=plotConfig$title, cex.main=1, cex.axis=0.8)
	#  	daysText <- paste("Day", i)
	#  	mtext(daysText, 3, cex=0.8)
	#  	# Plot offsetData points and actual offsetData lines
	#  	lines(offsetTimes, offsetData, col='steelblue', lty=1)
	#  	points(truncTimes, truncData, col='black', pch=16)
	#  	lines(fineTimes, allEvalFine$multiInf, lty=1)
	#  	# multiInfCoarse <- allEval$multiInf
	#  	multiInf <- allEvalFine$multiInf
	#  	for(k in 1:(length(allEvalFine$subInf))) {
	#  		sub <- allEvalFine$subInf[[k]]
	#  		subParams <- allEvalFine$subParams[[k]]
	#  		# Print sub epidemic graph
	#  		lines(fineTimes, sub, col=cl[k], lty=2)
	#  		lines(fineTimes, multiInf, col='black')
	#  		# lines(offsetTimes, multiInfCoarse, col='green')
	#  	}
	# 	legendText <- c("Epidemic Model", "Data")
	# 	legend("topright",legendText, col=c("black", "steelblue"), lty=1, cex=0.8)
	#  	dev.off()
	# }

	# # Set values of eval
	eval$multiParams <- optimParams
	eval$multiConds <- optimConds
	eval$optimTime <- optimTime
	eval$optimTimes <- c(ts[1:(k - 1)], optimTime)
	eval$k <- k
	eval$optimRSquare <- optimRSquare
	eval$epiTypes <- epiTypes
	# Get final residual from allEval infectious
	eval$finalResidual <- (offsetData[i] - allEval$multiInf[i])
	# Get all residuals up to current time
	eval$residuals <- (offsetData[1:i]) - (allEval$multiInf[1:i])
	# Set different multi and sub evals
	eval$nextPred <- allEval$multiInf[i + 1]
	eval$pastEval <- optimPastEval
	eval$allEval <- allEval
	eval$allEvalFine <- allEvalFine
	eval
}