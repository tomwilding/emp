orderOf <- function(x) {
	exp <- log(x, 10)

	exp <- floor(exp)
	10^exp
}