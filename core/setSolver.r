setSolver <- function(optimMethod, k) {
	# Select optimisation method
	parscale <- rep(c(-1,-1,1), k)
	switch(optimMethod,
		LMS = {
			optimSIRMulti <- function(times, data, initConds, initParams, ts, k) {
				params <- optim(initParams, sseSIRMulti, time=times, data=data, initConds=initConds, ts=ts, k=k, method="Nelder-Mead", control=list(parscale=parscale))
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