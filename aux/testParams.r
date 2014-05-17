testParams <- function(times, data, initConds, params, epiTypes, ts, k, granularity) {
	# I0 <- 1
	initConds <- c(100,1,0,0, 100,1,0,0)
	initParams <- c(-6.218659, -3.841935,  5.560658, 34.800554, -9.911200, -3.780002, 5.812445, 97.570066)
	epiTypes <- c(0, 4, 4)
	k <- 3
	granularity <- 1
	# times <- seq(from=0,to=60,by=0.02)
	eval <- evalMulti(times, data, initConds, initParams, epiTypes, k, granularity)
	plot(times, eval$multiInf, type="l")
	readline()
}