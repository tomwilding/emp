rSquareError <- function(predInfectious, data) {
	ssRes <- sum((predInfectious-data)^2)
	meanVec <- array(1,length(data))*myMean(data)
	ssTotal <- sum((meanVec-data)^2)
	if (ssTotal != 0) {
		rSqaure <- 1 - ssRes/ssTotal
	} else {
		rSquare <- 0
	}
}