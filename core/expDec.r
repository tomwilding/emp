expDec <- function(time, initConds, params) {
	i <- initConds[1]
	gamma <- exp(params[1])
	dI <- -gamma*i
	dT <- 0
	list(c(dI, dT))
}