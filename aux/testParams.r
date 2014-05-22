testParams <- function(times, data, initConds, params, epiTypes, ts, k, granularity) {
	# I0 <- 1
	initConds <- c(	1,1,0,0, 
				1,1,0,0)
	initParams <- c(-2.4010468, -1.3415407, -0.1093891,  0.4326234, -1.9578302, 1.3021430,  2.7858857, -0.7059237)
	epiTypes <- c(0, 4, 4)
	k <- 3
	ts <- c(1,13,52)
	# times <- seq(from=0,to=60,by=0.02)
	eval <- evalMulti(times, data, initConds, initParams, epiTypes, ts, k, granularity)
	fineTimes <- breakTime(times, granularity)
	plot(fineTimes, eval$multiInf, type="l")
	readline()
}