logit <- function(tmin, x, tmax) {
	if (x <= tmin) { 
		l <- tmin 
	} else if (x >= tmax) {
		l <- tmax
	} else {
		l <- log((x - tmin) / (tmax - x))
	}
	l
}