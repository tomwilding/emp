expDec <- function(time, initConds, params) {
	i <- initConds[1]
	gamma <- logisticTransform(1e-3, params[1], 1e-1)
	dI <- -gamma*i
	list(c(dI))
}