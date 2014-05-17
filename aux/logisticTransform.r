logisticTransform <- function(x, tmax) {
	tmax / (1 + exp(-x))
}