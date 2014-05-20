fitInRangeParallel <- function(optimSIRMulti, i, offsetTimes, offsetData, initConds, initParams, epiTypes, ts, k, plotConfig, p) {
	require(doMC)

	# Ensure length of range is greater than 0
	# assert_all_are_true(length(range) > 0)

	# Define eval vector
	fit <- c()

	# Register multi core backend with n cores
	registerDoMC(min(length(range), 30))

	# Truncate offsetData set
	truncTimes <- offsetTimes[1:i]
	truncData <- offsetData[1:i]

	# # SSE plot
	# points <- 100
	# betaVals <- log(1*10^-seq(4,2,length=points))
	# gammaVals <- log(1*10^-seq(2,0,length=points))
	# s0Vals <- rep(log(762),points)
	# # s0Vals <- log(seq(600,800,length=points))
	# sseB <- c(); sseS0 <- c();
	# sse <- matrix(1, points, points)
	# for (i in 1:length(betaVals)) {
	# 	# sseB[i] <- sseMulti(c(betaVals[i], log(1), log(500)), offsetTimes, offsetData, initConds, epiTypes, ts, k)
	# 	sseS0[i] <- sseMulti(c(log(0.00257), log(0.473), s0Vals[i]), offsetTimes, offsetData, initConds, epiTypes, ts, k)
	# 	for (j in 1:length(gammaVals)) {
	# 		sse[i, j] <- sseMulti(c(betaVals[i], gammaVals[j], s0Vals[i]), offsetTimes, offsetData, initConds, epiTypes, ts, k)
	# 	}
	# }
	# # plot(betaVals, sseB)
	# # plot(s0Vals, sseS0)
	# setEPS()
	# postscript(paste(plotConfig$fileName, "fluSurface.eps"), sep='')
	# persp(betaVals, gammaVals, sse, theta = -30, phi = 40, expand = 0.5, col = "lightblue", shade = 0.75, ticktype = "detailed",
 #     	xlab = "log(beta)", ylab = "log(gamma)", zlab = "SSE")
	# title(main="Optimisation surface over beta and gamma", cex.main=1, cex.axis=0.8)
	# dev.off()
	# print("ssePlotted")
	# readline()
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
	
	# time value t0 referenced from offset offsetData
	if (k > 1) {
		# Find optimal beta and gamma by optimising them to minimise the least square function
		# OptimSIRMulti passed in from call to setSolver
		tryCatch({
			optimParams <- optimSIRMulti(truncTimes, truncData, initConds, initParams, epiTypes, ts, k)
		}, warning = function(w) {
			print(w)
			print("optim warning")
		}, error = function(e) {
			print(e)
			print("optim failed")
		})
	}


	pastEval <- evalMulti(truncTimes, truncData, initConds, optimParams, epiTypes, ts, k, 1)
	predInfectiousPast <- pastEval$multiInf
	# rSquare error to determine best time to start fitting
	rSquareError <- rSquareError(predInfectiousPast, truncData)

	# TODO: Don't want to restrict eval - eval over all points in evalMulti after optim over all but last n
	allEvalFine <- evalMulti(offsetTimes, offsetData, initConds, optimParams, epiTypes, ts, k, timeStep)
	# Evaluate over all time
	allEval <- evalMulti(offsetTimes, offsetData, initConds, optimParams, epiTypes, ts, k, 1) 
	startOffset <- offsets$startOffset
	endOffset <- offsets$endOffset

	# Plot inline for dev
	if (p) {
	 	fineTimes <- breakTime(offsetTimes, timeStep)
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
	fit$initConds <- initConds
	# fit$optimTimes <- c(1, 23)
	fit$k <- k
	fit$optimRSquare <- rSquareError 
	# Get final residual from allEval infectious
	fit$finalResidual <- (offsetData[i] - allEval$multiInf[i])
	# Get all residuals up to current time
	fit$residuals <- (offsetData[1:i]) - (allEval$multiInf[1:i])
	# Set different multi and sub fits
	fit$nextPred <- allEval$multiInf[i+1]
	# fit$pastEval <- optimPastEval
	fit$allEval <- allEval
	fit$allEvalFine <- allEvalFine
	fit
}