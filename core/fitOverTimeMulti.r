fitOverTimeMulti <- function(optimMethod, times, data, initConds, initParams, epiTypes, offsets, thresholds, plotConfig) {
	
	# Unpack starting parameters, conditions and offsets
	startOffset <- offsets$startOffset
	endOffset <- offsets$endOffset

	startParams <- initParams
	startConds <- initConds
	# Take data set within specified offset
	offsetTimes <- times[startOffset:(length(times)-endOffset)]
	offsetData <- data[startOffset:(length(data)-endOffset)]

	# Max and min truncated data set sizes within offset data
	minTruncation <- offsets$minTruncation
	maxTruncation <- length(offsetData)
	# Min and max t0 values to explore within offset data
	minTRange <- 3
	maxTRange <- 3
	
	# Initialise other parameters
	rSquare <- 0
	rSquareRefit <- 0 
	totalRSqaure <- 0
	# Step size for iterative fitting
	step <- 1
	# Initial t0 value
	ts <- c(1)
	# Set the number of epidemics
	k <- 1

	# All evaluation vector
	evalList <- c()
	
	# Number of increasing residuals
	nResiduals <- 2
	residuals <- c()
	
	################################################# Decompose Epidemics ################################################
	# Truncate the data to i data points from minTruncation within offset data
	for (i in seq(from=minTruncation, to=maxTruncation, by=step)) {
		# Fit k epidemics
		print("### Fit k", quote=FALSE); print(paste(c("fitting "," of "), c(i, maxTruncation)), quote=FALSE)
		epidemicType <- epiTypes[k]
		if (epidemicType == 3) {
			eval <- fitInRangeParallel(setSolver(optimMethod, k, epiTypes), i, offsetTimes, offsetData, initConds, initParams, epiTypes, ts, k, c(minTRange:(i-maxTRange)), plotConfig)
		} else if (epidemicType == 1) {
			eval <- fitInRangeParallel(setSolver(optimMethod, k, epiTypes), i, offsetTimes, offsetData, initConds, initParams, epiTypes, ts, k, c(ts[k]:ts[k]), plotConfig)
		}
		# Store eval
		evalList[[i]] <- eval

		maxt <- eval$optimTime
		ts[k] <- maxt
		rSquare <- eval$optimRSquare
		multiParams <- eval$multiParams
		# Get last residual and update last n residuals vector
		residuals <- c(residuals, eval$finalResidual)


		# Try to improve fit if rSquare has deteriorated
		lim <- thresholds$lim
		diff <- thresholds$diff
		epidemicType <- getEpidemicType(residuals, nResiduals, rSquare)
		if ((rSquare < lim) && (epidemicType > 0)) {
			# Set new epidemic type in epidemic type array
			epiTypes <- c(epiTypes, epidemicType)
			# Try k+1 epidemics
			print(">>> Fit k+1", quote=FALSE)
			print(paste("!!!Epidemic", epidemicType), quote=FALSE)
			# Set initial parameters
			initParamsMulti <- newParams(initParams, startParams, epidemicType)
			initCondsMulti <- newParams(initConds, startConds, epidemicType)
			if (epidemicType == 3) {
				# Fit k+1 epidemics with new SIR sub epidemic
				eval <- fitInRangeParallel(setSolver(optimMethod, k+1, epiTypes), i, offsetTimes, offsetData, initCondsMulti, initParamsMulti, epiTypes, ts, k+1, c(minTRange:(i-maxTRange)), plotConfig)
			} else if (epidemicType == 1) {
				# Fit k+1 epidemics with set t0 at i
				eval <- fitInRangeParallel(setSolver(optimMethod, k+1, epiTypes), i, offsetTimes, offsetData, initCondsMulti, initParamsMulti, epiTypes, ts, k+1, c((i-1):(i-1)), plotConfig)
			}
			# TODO: Does the residuals array need to be updated here? = Only consider residuals in initial fitInRangeParallel above - will be updated in next loop
			# residuals <- nIncResiduals[(i-(nResiduals-1)):i]
			multiRSquare <- eval$optimRSquare
			# Set k+1 epidemics from now on
			k <- k + 1
			# Update parameters to continue fitting with k+1 epidemics
			multiParams <- eval$multiParams
			maxt <- eval$optimTime
			ts <- c(ts,maxt)
			rSquare <- multiRSquare
			initConds <- initCondsMulti
			# }
		}
		totalRSqaure <- totalRSqaure + rSquare
		# Update the initial parameters for the next fitting
		initParams <- multiParams
	}
	
	# Save all params
	save(evalList, file=plotConfig$dataFile)
	avRS <- totalRSqaure / length(seq(from=minTruncation, to=maxTruncation, by=step))
	print(avRS)
	avRS
	# save.image(plotConfig$envFile)
}

# Compare the last n residuals to sd of previous residuals before them
getEpidemicType <- function(residuals, nRes, rSquare) {
	incRes <- c()
	type <- 0
	# Standard deviation of residuals
	sdRes <- 0
	# Get the last n residuals
	residuals <- (residuals[!is.na(residuals)])
	# If less then one residual then set residual check to false
	minIncRes <- Inf
	
	# Ensure more than one residual before the last n residuals to calculate sdRes
	resLength <- length(residuals)
	print
	if (resLength > nRes + 1) {
		# Get standard deviation of residuals before the ones considered
		meanRes <- myMean(residuals[(1:(resLength - nRes))])
		sdRes <- mySd(residuals[1:(resLength - nRes)], meanRes)
		# Reset between epidemics
		meanDiffRes <- myMeanDiff(residuals[1:(resLength - nRes)])
		sdDiffRes <- mySdDiff(residuals[1:(resLength - nRes)], meanDiffRes)
		diffRes <- abs(residuals[resLength - nRes + 1] - residuals[resLength - nRes])
		print(paste("sdRes ", sdRes))
		print(paste("RSq", rSquare))
		print(paste("MeanDiffRes", meanDiffRes))
		print(paste("SdDiffRes", sdDiffRes))
		print(paste("DiffRes", diffRes))
		# If current residual sd is above zero check if last n residuals are above set number of sd
		if (sdRes > 0) {
			# Index of first residual to check
			startResIndex <- resLength - nRes
			# Assume incRes is True and check condition for all n residuals
			# sameSign <- max(residuals[(startResIndex + 1):resLength]) < 0 || min(residuals[(startResIndex + 1):resLength]) > 0
			sigIncRes <- diffRes > (meanDiffRes + 2 * sdDiffRes)

			incResiduals <- TRUE
			for (i in 1:nRes) {
				# Check magnitude of residuals
				print(residuals[startResIndex + i])
				minIncRes <- min(minIncRes, (abs(residuals[startResIndex + i])))
			}
		}
		# print(paste("SameSign",sameSign))
		# Set epidemic type according to residual limit
		sirSD <- sdRes * 8
		spikeSD <- sdRes * 2
		# If minimum residual increase is more than required, then set type
		if (minIncRes > spikeSD && sigIncRes) { #&& sameSign) {
			print("T1 Set")
			# readline()
			type <- 1
		} else if (minIncRes > sirSD && sigIncRes) { #&& sameSign) {
			print("T2 Set")
			# readline()
			type <- 3
		}
	}
	type
}

newParams <- function(initVec, startVec, epidemicType) {
	if (epidemicType == 1) {
		# Update params for spike epidemic take I0 and gamma from initial
		newVec <- c(initVec, startVec[2])
	} else if (epidemicType == 3) {
		# Update params for SIR epidemic
		newVec <- c(initVec, startVec)
	}

}