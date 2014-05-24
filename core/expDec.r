expDec <- function(time, initConds, params) {
	i <- initConds[1]
	gamma <- logisticTransform(1e-3, params[1], 0.5)
	dI <- -gamma*i
	list(c(dI))
}