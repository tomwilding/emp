logit <- function(x, maxt) {
	if (x == 1) {
		l <- 1
	} else {
		l <- log(x / (1 - maxt*x))
	}
}