logistic <- function(tmin, x, tmax) {
	((tmax - tmin) / (1 + exp(-x))) + tmin
}