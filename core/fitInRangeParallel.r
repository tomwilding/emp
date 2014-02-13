fitInRangeParallel <- function(optimSIRMulti, i, times, data, initConds, initParams, ts, k, range, actualFit, plotConfig) {
	require(doMC)
	# source('breakTime.r')

	# cl <- makeCluster(3)
	# registerDoParallel(cl)
	# clusterExport(cl, c("sseSIRMulti", "lsoda"))
	# numCores <- 15
	# registerDoMC(numCores)
	# print(length(range))
	registerDoMC(length(range))
	# Truncate data set
	truncTimes <- times[1:i]
	truncData <- data[1:i]

	# Fine Times for evaluation
	timeStep <- 0.05
	fineTimes <- breakTime(times, timeStep);

	rSquareMulti <- 0
	maxRSquare <- 0
	optimParams <- initParams

	# Range of feasible t0 values to explore referenced from offset data
	# t0 <- proc.time()
	allOptimParams <- foreach (t=range, .combine="rbind") %dopar% {
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

		predInfectiousAll <- evalSIRMulti(times, data, initConds, optimParams, tsExplore, k)
		predInfectiousPast <- predInfectiousAll[1:i]
		predInfectiousFuture <- predInfectiousAll[(i+1):length(predInfectiousAll)]
		if ((i+8) < length(predInfectiousAll)) {
			predInfectiousFuture1 <- predInfectiousAll[(i+1):(i+8)]
			futureData1 <- data[(i+1):(i+8)]
		} else {
			futureData1 <-  data[(i+1):length(data)]
		}
		# if (i==8) {
		# 	print(futureData1)
		# 	print(predInfectiousFuture1)
		# }
		# TODO: Don't calc all Rsq inside parallel section
		futureData <- data[(i+1):length(data)]
		rSquareErrorAll <- rSquareError(predInfectiousAll, data)
		rSquareErrorPast <- rSquareError(predInfectiousPast, truncData)
		rSquareErrorFuture <- rSquareError(predInfectiousFuture, futureData)
		rSquareErrorFuture1 <- rSquareError(predInfectiousFuture, futureData1)
		finalRes <- (data[i] - predInfectiousAll[i])
		c(rSquareErrorFuture, rSquareErrorFuture1, rSquareErrorAll, finalRes, t, rSquareErrorPast, optimParams)
	}
	# t1 <- proc.time()
	# print(t1-t0)
	# readline()
	# Find max rSqaureError entry
	if (length(range) > 1) {
		optimIndex <- which.max(allOptimParams[,6])
		optimParams <- allOptimParams[optimIndex,]
	} else {
		optimParams <- allOptimParams
	}
	rSquareFuture <- optimParams[1]
	rSquareFuture1 <- optimParams[2]
	rSquareAll <- optimParams[3]
	finalRes <- optimParams[4]
	rSquarePast <- optimParams[6]
	optimParams <- optimParams[5:length(optimParams)]

	# EvaluateSIR over to plot fitted epidemics
	par(mar=c(7.1,4.1,4.1,2.1))
	plot(1:length(times), data, xlab='Epochs', ylab='Infected Individuals', col='steelblue')
	title(main=plotConfig$title, cex.main=1, cex.axis=0.8)
	daysText <- paste("Epochs after outbreak = ", i)
	mtext(daysText, 3, cex=0.8)
	RSqPastText <- paste("Past RSquare = ", signif(rSquarePast, digits=3))
	RSqFutureText <- paste("Future RSquare = ", signif(rSquareFuture, digits=3))
	# RSqFuture1Text <- paste("Future T7 RSquare = ", signif(rSquareFuture1, digits=3))
	RSqAllText <- paste("All RSquare = ", signif(rSquareAll, digits=3))
	# mtext(RSqText, 1, at=10, padj=6, cex=0.5)
	mtext(RSqPastText, 1, at=plotConfig$rat, padj=6, cex=0.7) #other
	mtext(RSqFutureText, 1, at=plotConfig$rat, padj=8, cex=0.7) #other
	# mtext(RSqFuture1Text, 1, at=plotConfig$rat, padj=10, cex=0.7) #other
	mtext(RSqAllText, 1, at=plotConfig$rat, padj=10, cex=0.7) #other
	# mtext(RSqText, 1, at=30, padj=6, cex=0.5) #bl
	# mtext(RSqText, 1, at=20, padj=6, cex=0.5) #hs
	lines(1:length(times), data, col='steelblue')
	points(1:length(truncTimes), truncData, col='black', pch=16)
	allPredInf <- decomposeEpidemics(times, data, initConds, optimParams[3:length(optimParams)], c(ts[1:k-1],optimParams[1]), k, actualFit, plotConfig)
	allPredInf$optimParams <- optimParams
	allPredInf$rSquares <- c(rSquareFuture, rSquareFuture1, rSquareAll, rSquarePast)
	allPredInf$finalRes <- finalRes;
	allPredInf
}