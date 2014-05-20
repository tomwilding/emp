expDec <- function(time, initConds, params) {
	i <- initConds[1]
	gamma <- logisticTransform(1e-2, params[1], 0.1)
	dI <- -gamma*i
	list(c(dI))
}