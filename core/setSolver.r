setSolver <- function(optimMethod, k, epiTypes) {
	# require(neldermead)
	# print("Set Solver")
	# Set parscale for optimisation
	# parscale <- c()
	# for (t in epiTypes) {
	# 	if (t == 1) {
	# 		parscale <- c(parscale, c(1))
	# 	} else if (t == 4) {
	# 		parscale <- c(parscale, c(-1,-1,1,1))
	# 	}
	# }

	# Select optimisation method
	switch(optimMethod,
		LMS = {
			optimSIRMulti <- function(times, data, initConds, initParams, epiTypes, ts, k, timeStep) {
				optimisationParameters <- initParams
				for (o in 1 : 10) {
					print(paste("optim",o))
					params <- optim(optimisationParameters, sseMulti, time=times, data=data, initConds=initConds, epiTypes=epiTypes, ts=ts, k=k, timeStep=timeStep, control=list(maxit=1000))
					optimisationParameters <- params$par
					print(params)
				}
				# if (params$convergance > 0) {

				# }
				optimisationParameters
			}
		},
		SP = {
			optimSIRMulti <- function(times, data, initConds, initParams, epiTypes, ts, k, timeStep) {
				params <- optim(initParams, sseMulti, time=times, data=data, initConds=initConds, epiTypes=epiTypes, ts=ts, k=k, timeStep=timeStep, control=list(maxit=1000))
				optimisationParameters <- params$par
			}
		},
		MLE = {
			optimSIRMulti <- function(times, data, initConds, initParams, epiTypes, ts, k) { 
				namedParams <- nameParams(initParams)
				mle <- mle2(sirNegLL, start=namedParams, data=list(timeIn=times,  dataIn=data, initConds=initConds, ts=ts, k=k, epiTypes=epiTypes), method="Nelder-Mead")
				params <- as.list(coef(mle))
				optimParams <- unnameParams(params)
				# params <- optim(initParams, sirNegLL, time=times, data=data, initConds=initConds, ts=ts, k=k, epiTypes=epiTypes, method="Nelder-Mead", control=list(parscale=parscale))
				# optimParams <- params$par
			}
		},
		{
			print('Solver method not specified or not recognised, defaulting to LMS')
			optimSIRMulti <- function(times, data, initConds, initParams, ts, k) { 
				params <- optim(initParams, sseMulti, time=times, data=data, initConds=initConds, ts=ts, k=k, epiTypes=epiTypes, method="Nelder-Mead", control=list(parscale=parscale))
				optimParams <- params$par
			}
		}
	)
}

# nameParams <- function(initParams) {
# 	list(b=initParams[1], g=initParams[2], s0=initParams[3])
# }

# unnameParams <- function(params) {
# 	c(params$b, params$g, params$s0)
# }

# sirNegLL <- function(b , g, s0, timeIn, dataIn, initConds, ts, k, epiTypes) {
# 	granularity <- 1
# 	params <- c(b,g,s0)
# 	gamma <- exp(g)
# 	if (gamma > 1 || gamma <= 1e-4) {
# 		nll <- -Inf
# 	} else {
# 		eval <- evalMulti(timeIn, dataIn, initConds, params, epiTypes, ts, k, granularity)
# 		nll <- -sum(dnorm(x=dataIn, mean=eval$multiInf, sd=1, log=TRUE))
# 	}
# 	nll
# }
# sirNegLL <- function(params, times, data, initConds, epiTypes, ts, k) {
# 	granularity <- 1
# 	eval <- evalMulti(times, data, initConds, params, epiTypes, ts, k, granularity)
# 	nll <- -sum(dpois(x=data, lambda=eval$multiInf, log=TRUE))	
# }

# myOptim <- function(initParams, sseMulti, times, data, initConds, ts, k) {
# 	# Create parscale array for each sub epidemic parameter set
# 	parscale <- rep(c(-1,-1,1), k)
# 	# Flatten lists of parameters to send to optim
# 	initParamsFlatten <- flatten(initParams)
# 	initCondsFlatten <- flatten(initConds)
# 	# Optimise over parameter array
# 	params <- optim(initParamsFlatten$flattenList, sseMulti, time=times, data=data, initConds=initCondsFlatten$flattenList, ts=ts, k=k, method="Nelder-Mead", control=list(parscale=parscale))
# 	print(unflatten(params$par))
# 	optimParams <- unflatten(params$par)
# }

# flatten <- function(list) {
# 	flatten <- c()
# 	flattenList <- c()
# 	subSize <- c()
# 	# Loop through all lists
# 	for (i in 1:length(list)) {
# 		innerList <- list[[i]]
# 		# Flatten list to array
# 		flattenList <- c(flattenList, innerList)
# 		# Record size of sub epidemics parameters in array of sizes
# 		subSize[i] <- length(innerList)
# 	}
# 	flatten$subSize <- subSize
# 	flatten$flattenList <- flattenList
# 	flatten
# }

# unflatten <- function(list, subSize) {
# 	newList <- list()
# 	newList[[1]] <- list[1:3]
# 	newList
# }