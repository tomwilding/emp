testParams <- function(times, data, initConds, params, epiTypes, ts, k, granularity) {
	# I0 <- 1
	initConds <- c(	1,1,0,0, 
				1,1,0,0)
	initParams <- c(32.79236020,  23.66425693, -60.66457472, 4.51213643,  -2.01868608,
					-26.15104230,  -0.02743762,  34.40602921)
	epiTypes <- c(0, 4, 4)
	k <- 3
	granularity <- 1
	# times <- seq(from=0,to=60,by=0.02)
	eval <- evalMulti(times, data, initConds, initParams, epiTypes, ts, k, granularity)
	plot(times, eval$multiInf, type="l")
	readline()
}