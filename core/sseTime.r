# Helper function for LMS to compute sse
sseTime <- function(startTimes, times, data, initConds, initParams, epiTypes, ts, k, timeStep) {
	print(startTimes)
	params <- optim(initParams$params, sseMulti, time=times, data=data, initConds=initConds, epiTypes=epiTypes, ts=ts, startTimes=startTimes, k=k, timeStep=timeStep)
	optimisationParameters <- params$par
	print(optimisationParameters)
	initParams$params <- optimisationParameters
	eval <- evalMulti(times, data, initConds, optimisationParameters, epiTypes, ts, startTimes, k, timeStep)
	# print(eval$multiParams)
	predInf <- getObservations(eval$multiInf, timeStep)
	# eval <- evalMulti(times, data, initConds, params, epiTypes, ts, k, 1)
	# predInf <- eval$multiInf
	sse <- sum((predInf - data)^2)
}