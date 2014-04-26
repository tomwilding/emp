rSquareError <- function(predInfectious, data) {
	ssRes <- sum((predInfectious-data)^2)
	meanVec <- array(1,length(data))*myMean(data)
	ssTotal <- sum((meanVec-data)^2)
	# ssReg <- sum((meanVec-predInfectious)^2)
	# rSquare1 <- ssReg/ssTotal
	rSquare <- 1 - ssRes/ssTotal
	# print(rSquare)
	# print(rSquare1)
}