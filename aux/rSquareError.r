rSquareError <- function(predInfectious, data) {
	ssRes <- sum((predInfectious-data)^2)
	meanVec <- matrix(1,length(data),1)*myMean(data)
	ssTotal <- sum((meanVec-data)^2)
	rSqaure <- 1 - ssRes/ssTotal
}