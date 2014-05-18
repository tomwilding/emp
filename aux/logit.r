logit <- function(tmin, x, tmax) {
	l <- log((x - tmin) / (tmax - x))
}