ssError <- function(predInfectious, data) {
	ssRes <- sum((predInfectious-data)^2)
}