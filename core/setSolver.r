setSolver <- function(optimMethod, k, epiTypes) {
	# require(bbmle)
	# print("Set Solver")
	# Set parscale for optimisation
	parscale <- c()
	for (t in epiTypes) {
		if (t == 1) {
			parscale <- c(parscale, c(1))
		} else if (t == 3) {
			parscale <- c(parscale, c(-1,-1,1))
		}
	}

	# Select optimisation method
	switch(optimMethod,
		LMS = {
			optimSIRMulti <- function(times, data, initConds, initParams, epiTypes, ts, k) {
				optimisationParameters <- initParams
				for (i in 1 : 5) {
					params <- optim(optimisationParameters, sseMulti, time=times, data=data, initConds=initConds, ts=ts, k=k, epiTypes=epiTypes)
					# print(params)
					optimisationParameters <- params$par
				}
				optimisationParameters
			}
		},
		MLE = {
			optimSIRMulti <- function(times, data, initConds, initParams, epiTypes, ts, k) { 
				namedParams <- nameParams(initParams)
				mle <- mle2(sirNegLL, start=namedParams, data=list(times=times,  data=data, initConds=initConds, ts=ts, k=k, epiTypes=epiTypes), method="Nelder-Mead", control=list(maxit=5000))
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

nameParams <- function(initParams) {
	list(b=initParams[1], g=initParams[2], s0=initParams[3])
}

unnameParams <- function(params) {
	c(params$b, params$g, params$s0)
}

sirNegLL <- function(b, g, s0, times, data, initConds, epiTypes, ts, k) {
	params <- c(b, g, s0)
	granularity <- 1
	eval <- evalMulti(times, data, initConds, params, epiTypes, ts, k, granularity)
	nll <- -sum(dpois(x=data, lambda=eval$multiInf, log=TRUE))	
}

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