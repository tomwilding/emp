testParams <- function(times, data, initConds, params, epiTypes, ts, k, granularity) {
	# I0 <- 1
	initConds <- c(1000,1,0,0, 1000,1,0,0)
	initParams <- c(-10.8340994,  -3.8988462,   8.8141918 ,0.3772866,  -7.6489684,  -4.6050914, 6.9343062, -1.0678850)
	epiTypes <- c(0, 4, 4)
	k <- 3
	granularity <- 1
	# times <- seq(from=0,to=60,by=0.02)
	eval <- evalMulti(times, data, initConds, initParams, epiTypes, ts, k, granularity)
	plot(times, eval$multiInf, type="l")
	readline()
}