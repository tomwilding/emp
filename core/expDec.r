expDec <- function(time, initConds, params) {
	i <- initConds[1]
	gamma <- exp(params[1])
	dI <- -gamma*i
	list(c(dI))
}