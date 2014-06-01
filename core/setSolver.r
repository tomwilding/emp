setSolver <- function(optimMethod, k, epiTypes) {

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
				optimisationParameters
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