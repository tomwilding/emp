fitOverTimeMulti <- function(optimMethod, times, data, initConds, initParams, epiTypes, offsets, target, plotConfig) {

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
	k <- 1
	
	# Epidemic start times
	ts <- c()

	# Eval store	 
	evalList <- c()

	################################################# Decompose Epidemics ################################################
	# Truncate the data to i data points from minTruncation within offset data
	for (i in seq(from=minTruncation, to=maxTruncation, by=1)) {
		# gradientSearch(times[1:i], data[1:i])
		
		print("### Fit k", quote=FALSE); print(paste(c("fitting "," of "), c(i, maxTruncation)), quote=FALSE)
		if ( k > 2 ) {
			write(paste(c("fitting "," of "), c(i, maxTruncation)), file="optimParams.txt", append=TRUE)
			write("", file="optimParams.txt", append=TRUE)
		}

		# Optimise k epidemics
		optimParams <- initParams
		for (rep in 1 : 10) {
			eval <- fitInRangeParallel(setSolver(optimMethod, k, epiTypes), i, offsetTimes, offsetData, initConds, optimParams, epiTypes, k, ts)
			optimParams <- eval$multiParams
			if ( k > 2 ) {write("optimLoop", file="optimParams.txt", append=TRUE)}
		}
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
		if ((rSquare > 0 && rSquare < target) || outbreakDetected) {
			print("Epi Detected")
			# readline()
			k <- k + 1
			# lastestEpidemicType <- determineEpidemicType()
			lastestEpidemicType <- 5;
			epiTypes[k] <- lastestEpidemicType
			# Update start times
			ts <- c(ts, i)
			# Update parameters
			initParams <- newParams(initParams, i, lastestEpidemicType, nResiduals, eval)
			initConds <- newConds(initConds, i, lastestEpidemicType, nResiduals)
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
		print(paste("lim", meanRes + 3*sdRes))
		# Check n residuals are outside of range
		for (i in 0:(nResiduals - 1)) {
			res <- residuals[length(residuals) - i]
			outbreak <- outbreak && (res > (meanRes + 3*sdRes))
		}
	}
	outbreak
}

# Compare the last n residuals to sd of previous residuals before them
determineEpidemicType <- function(residuals, nResisuals, window, rSquare) {
	incRes <- c()
	type <- 0
	# Standard deviation of residuals
	sdRes <- 0
	# Get the last n residuals
	residuals <- (residuals[!is.na(residuals)])
}

newParams <- function(initVec, i, epidemicType, nResisuals, eval) {
	initTime <- i - nResisuals
	initialInf <- data[initTime + 1]

	# Update parameters
	if (length(initVec > 0)) {
		initParams <- initVec[1:5]
		# Update start time
		initParams[5] <- i
		initParams <- c(initVec, initParams)
	# Initial parameters
	} else {
		initParams <- c(log(0.001), log(0.1), log(orderRange(initialInf)*10), log(orderRange(initialInf)), logit(initTime, i))
	}
	initParams
}

newConds <- function(initVec, i, epidemicType, nResisuals) {
	if (length(initVec > 0)) {
		# Update conds
		initConds <- initVec[1:5]
		initConds <- c(initVec, initConds)
	} else {
		initConds <- c(1, 1, 0, 0, 0)
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