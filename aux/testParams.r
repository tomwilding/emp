testParams <- function(times, data, initConds, params, epiTypes, ts, k, granularity) {
	# I0 <- 1
initConds <- c(	1,1,0,0, 
				1,1,0,0,
				1,
				1,1,0,0,
				1,1,0,0)
	initParams <- c(-2.5814881 , -3.7259482  , 7.8348258 ,  1.0509491 , 
					-2.0631394 , -2.5093139, 8.2421060 , 22.1483512 , 
					-1.6689080 , 
					-2.2164175 ,-16.5738898,   8.6088228, 35.7868721 ,  
					1.5066700,  -6.5856233  , 8.2203729  , 0.2537642)
	epiTypes <- c(0, 4, 4, 1, 4 ,4)
	k <- 6
	granularity <- 1
	# times <- seq(from=0,to=60,by=0.02)
	eval <- evalMulti(times, data, initConds, initParams, epiTypes, ts, k, granularity)
	plot(times, eval$multiInf, type="l")
	readline()
}