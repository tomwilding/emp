setSolver <- function(optimMethod, k) {
	# print("Set Solver")
	# Select optimisation method
	parscale <- rep(c(-1,-1,1), k)
	switch(optimMethod,
		LMS = {
			optimSIRMulti <- function(times, data, initConds, initParams, epiTypes, ts, k) {
				params <- optim(initParams, sseSIRMulti, time=times, data=data, initConds=initConds, ts=ts, k=k, epiTypes=epiTypes, method="Nelder-Mead", control=list(parscale))
				# myOptim(initParams, sseSIRMulti, times, data, initConds, ts, k)
				optimParams <- params$par
			}
		},
		MLE = {
			optimSIRMulti <- function(times, data, initConds, initParams, ts, k) { 
				namedParams <- nameParams(initParams)
				mle <- mle2(minuslogl=buildSirNegLLMulti(initParams), start=namedParams, data=list(timeIn=times,  dataIn=data, initConds=initConds, ts=ts, k=k), method="Nelder-Mead", control=list(parscale=c(-1,-1,1,-1,-1,1)))
				params <- as.list(coef(mle))
				optimParams <- unnameParams(params)
			}
		},
		{
			print('Solver method not specified or not recognised, defaulting to LMS')
			optimSIRMulti <- function(times, data, initConds, initParams, ts, k) { 
				params <- optim(initParams, sseSIRMulti, time=times, data=data, initConds=initConds, ts=ts, k=k, method="Nelder-Mead", control=list(parscale=c(-1,-1,1,-1,-1,1)))
				optimParams <- params$par
			}
		}
	)
}

# myOptim <- function(initParams, sseSIRMulti, times, data, initConds, ts, k) {
# 	# Create parscale array for each sub epidemic parameter set
# 	parscale <- rep(c(-1,-1,1), k)
# 	# Flatten lists of parameters to send to optim
# 	initParamsFlatten <- flatten(initParams)
# 	initCondsFlatten <- flatten(initConds)
# 	# Optimise over parameter array
# 	params <- optim(initParamsFlatten$flattenList, sseSIRMulti, time=times, data=data, initConds=initCondsFlatten$flattenList, ts=ts, k=k, method="Nelder-Mead", control=list(parscale=parscale))
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