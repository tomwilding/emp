expDec <- function(time, initConds, params) {
	i <- initConds[1]
	gamma <- exp(params[1])
	dI <- -gamma*i
	dt <- 0
	list(c(dI, dt))
}