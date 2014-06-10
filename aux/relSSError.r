relSSError <- function(predInfectious, data) {
	# print(abs(predInfectious-data) / data)
	absDiff <- abs(predInfectious - data)
	absDiff[data == 0] <- 0
	data[data == 0] <- 1
	ssRes <- sum(absDiff/ data)
}