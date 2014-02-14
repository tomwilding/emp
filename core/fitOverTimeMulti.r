fitOverTimeMulti <- function(optimMethod, times, data, initConds, initParams, offsets, thresholds, plotConfig) {
	# library(bbmle, lib.loc="Rpackages")
	# source('breakTime.r')
	# source('sseSIRMulti.r')
	# source('evalSIRMulti.r')
	# source('rSquareError.r')
	# source('fitInRangeParallel.r')
	# source('decomposeEpidemics.r')
	# source('setSolver.r')

	# Unpack starting parameters, conditions and offsets
	startOffset <- offsets$startOffset
	endOffset <- offsets$endOffset

	startParams <- initParams
	startConds <- initConds
	# Take data set within specified offset
	offsetTimes <- times[startOffset:(length(times)-endOffset)]
	offsetData <- data[startOffset:(length(data)-endOffset)]
	# Max and min truncated data set sizes within offset data
	maxTruncation <- length(offsetData)
	minTruncation <- offsets$minTruncation
	# Min and max t0 values to explore within offset data
	minTRange <- 3
	maxTRange <- 3
	
	# Initialise other parameters
	rSquare <- 0
	rSquareRefit <- 0 
	# Step size for iterative fitting
	step <- 1
	# Initial t0 value
	ts <- c(1)
	# Set the number of epidemics
	k <- 1
	allPredInfList <- c()
	# Number of increasing residuals
	nRes <- 3
	nIncRes <- c()

	################################################# Decompose Epidemics ################################################
	# Truncate the data to i data points from minTruncation within offset data
	for (i in seq(from=minTruncation, to=maxTruncation, by=step)) {
		# Set graph printing settings
		setEPS()
		graphName <- paste("t", i, sep='')
		graphName <- paste(graphName, ".eps", sep='')
		postscript(paste(plotConfig$fileName, graphName, sep=''))		

		# Fit k epidemics
		print("### Fit k", quote=FALSE); print(paste(c("fitting "," of "), c(i, maxTruncation)), quote=FALSE)
		allPredInf <- fitInRangeParallel(setSolver(optimMethod, k), i, offsetTimes, offsetData, initConds, initParams, ts, k, c(minTRange:(i-maxTRange)), 1, plotConfig)
		# Store allPredInf
		allPredInfList[[i]] <- allPredInf

		maxt <- allPredInf$optimTime
		ts[k] <- maxt
		rSquare <- allPredInf$optimRSquare
		optimParams <- allPredInf$optimParams


		# Try to improve fit if rSquare has deteriorated
		lim <- thresholds$lim
		diff <- thresholds$diff
		# if (rSquare < lim && incMag(lastNRes, nRes)) {
		if (rSquare < lim) {
			# Try k+1 epidemics
			print(">>> Fit k+1", quote=FALSE)
			# Set initial parameters
			initParamsMulti <- c(initParams, startParams)
			initCondsMulti <- c(initConds, startConds)
			# Fit k+1 epidemics
			allPredInf <- fitInRangeParallel(setSolver(optimMethod, k+1), i, offsetTimes, offsetData, initCondsMulti, initParamsMulti, ts, k+1, c(minTRange:(i-maxTRange)), 2, plotConfig)
			# nIncRes[i] <- allPredInf$finalRes
			# lastNRes <- nIncRes[(i-(nRes-1)):i]
			multiRSquare <- allPredInf$optimRSquare

			# If k+1 is significantly better then continue with k+1 fit
			# if ((multiRSquare - rSquare) > diff || incMag(lastNRes, nRes)) {
			if ((multiRSquare - rSquare) > diff) {
				print("!!!Set k+1", quote=FALSE)
				# Set k+1 epidemics from now on
				k <- k + 1
				# Update parameters to continue fitting with k+1 epidemics
				optimParams <- allPredInf$optimParams
				maxt <- allPredInf$optimTime
				ts <- c(ts,maxt)
				rSquare <- multiRSquare
				initConds <- initCondsMulti
			}
		}
		# Update the initial parameters for the next fitting
		initParams <- optimParams
	}
	# Save all params
	save(allPredInfList, file=plotConfig$dataFile)
	save.image(plotConfig$envFile)
}

incMag <- function(lastNRes, n) {
	incMag <- FALSE
	inc <- TRUE
	ssign <- TRUE
	lastNRes <- (lastNRes[!is.na(lastNRes)])
	if (length(lastNRes) == n) {
		# Check is continuous increasing magnitude
		for (i in 2:n) {
			inc <- inc && (abs(lastNRes[i]) > abs(lastNRes[i-1]))
			ssign <- ssign && ((lastNRes[i]<0) == (lastNRes[i-1]<0))
		}
		incMag <- inc && ssign
		print("IM") 
		print(incMag)
	}
	incMag
}