rSquareError <- function(predInfectious, data) {
	ssRes <- sum((predInfectious-data)^2)
	mean <- mean(data)
	ssTotal <- sum((mean-data)^2)
	rSquare <- 1 - ssRes/ssTotal
}