fitOverTimeMulti <- function(optimMethod, times, data, initConds, initParams, epiTypes, offsets, target, plotConfig, tmax) {
	
	# Unpack starting parameters, conditions and offsets
	startOffset <- offsets$startOffset
	endOffset <- offsets$endOffset

	# Take data set within specified offset
	offsetTimes <- times[startOffset:(length(times)-endOffset)]
	offsetData <- data[startOffset:(length(data)-endOffset)]

	# Max and min truncated data set sizes within offset data
	minTruncation <- offsets$minTruncation
	maxTruncation <- length(offsetData)
	
	# Initial t0 value
	# ts <- c(1, 23)
	# Set the number of epidemics
	k <- 3
	
	evalList <- c()

	################################################# Decompose Epidemics ################################################
	# Truncate the data to i data points from minTruncation within offset data
	for (i in seq(from=minTruncation, to=maxTruncation, by=1)) {
		print("### Fit k", quote=FALSE); print(paste(c("fitting "," of "), c(i, maxTruncation)), quote=FALSE)
		if ( k > 2 ) {
			write(paste(c("fitting "," of "), c(i, maxTruncation)), file="optimParams.txt", append=TRUE)
			write("", file="optimParams.txt", append=TRUE)
		}

		# Optimise k epidemics
		optimParams <- initParams
		for (rep in 1 : 10) {
			eval <- fitInRangeParallel(setSolver(optimMethod, k, epiTypes), i, offsetTimes, offsetData, initConds, optimParams, epiTypes, k, tmax)
			optimParams <- eval$multiParams
			if ( k > 2 ) {write("optimLoop", file="optimParams.txt", append=TRUE)}
		}
		# maxt <- eval$optimTime
		# ts[k] <- maxt
		rSquare <- eval$optimRSquare
		multiParams <- eval$multiParams
		print(multiParams)
		# print(paste("R0", (exp(multiParams[1])*exp(multiParams[3])) / exp(multiParams[2])))
		print(paste("rs",rSquare))
		# TODO: Store eval starts at i index causing NA
		evalList[[i]] <- eval

		# Check for redundant epidemics

		# # If rSquare has deteriorated then fit k + 1 epidemics
		# if (rSquare < target) {
		# 	# Update k
		# 	k <- k + 1
			
		# 	# Set parameters

		# Optimise k + 1 epidemics
		nResiduals <- 2
		outbreakDetected <- detectOutbreak(eval$residuals, nResiduals)
		if ((rSquare > 0 && rSquare < target || outbreakDetected)) {
			print("Epi Detected")
			readline()
			k <- k + 1
			lastestEpidemicType <- 4;
			epiTypes[k] <- lastestEpidemicType
			# Update parameters
			print(paste("ip",initParams))
			print(paste("ic", initConds))
			initParams <- newParams(initParams, i, lastestEpidemicType)
			initConds <- newConds(initConds, i, lastestEpidemicType)
		}

		# }
		# Set start parameters to current optimised parameters
		# initParams <- multiParams
	}
	
	# Save all params
	save(evalList, file=plotConfig$dataFile)
	# save.image(plotConfig$envFile)
}

detectOutbreak <- function(residuals, nResiduals) {
	outbreak <- FALSE
	if (length(residuals) > nResiduals) {
		outbreak <- TRUE
		pastResiduals <- residuals[1 : (length(residuals) - nResiduals)]
		meanRes <- myMean(pastResiduals)
		sdRes <- mySd(pastResiduals, meanRes)
		# print(meanRes)
		# print(sdRes)
		print(meanRes + 3*sdRes)
		# Check n residuals are outside of range
		for (i in 0:(nResiduals - 1)) {
			print(residuals[length(residuals) - i])
			outbreak <- outbreak && (residuals[length(residuals) - i] > (meanRes + 3*sdRes))
		}
	}
	outbreak
}

# Compare the last n residuals to sd of previous residuals before them
getEpidemicType <- function(residuals, nRes, window, rSquare) {
	incRes <- c()
	type <- 0
	# Standard deviation of residuals
	sdRes <- 0
	# Get the last n residuals
	residuals <- (residuals[!is.na(residuals)])


}

newParams <- function(initVec, i, epidemicType) {
	if (length(initVec > 0)) {
		# Update params for SIR epidemic
		initParams <- initVec
		initParams[4] <- logit(i, tmax)
		initParams <- c(initVec, initParams)
	} else {
		initParams <- c(log(0.001), log(0.1), log(data[i]*10), log(data[i]), logit(i, tmax))
	}
	print(initParams)
	initParams
}

newConds <- function(initVec, i, epidemicType) {
	if (length(initVec > 0)) {
		# Update conds
		initConds <- initVec
		initConds <- c(initVec, initConds)
	} else {
		initConds <- c(1, data[i], 0, 0, 0)
	}
	initConds
}

reduceParams <- function(initVec, epidemicType) {
	if (epidemicType == 1) {
		# Update params for spike epidemic used I0 and gamma from initial
		newVec <- initVec[1:(length(initVec) - 1)]
	} else if (epidemicType == 3) {
		# Update params for SIR epidemic
		newVec <- initVec[1:(length(initVec) - 3)]
	}
}

countStartTime <- function(startTime, startTimePrev, startTimeCount) {
	if (startTime == startTimePrev) {
		startTimeCount <- startTimeCount + 1
	}
	startTimeCount
}