testParams <- function(times, data, initConds, params, epiTypes, ts, k, granularity) {
	# I0 <- 1
	initConds <- c(1,1,0)
	initParams <- c(-6.382314, -4.605170,  5.578618)
	epiTypes <- c(0, 3)
	ts <- c(1, 30)
	k <- 2
	granularity <- 1
	# times <- seq(from=0,to=60,by=0.02)
	eval <- evalMulti(times, data, initConds, initParams, epiTypes, ts, k, granularity)
	plot(times, eval$multiInf, type="l")
	readline()
}