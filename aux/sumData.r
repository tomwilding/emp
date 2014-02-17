sumData <- function(data, n) {
	takeData <- c()
	for (i in 1:(length(data)%/%n)) {
		startIndex <- ((((i-1)*n)+1))
		endIndex <- (i*n)
		binData <- data[startIndex:endIndex]
		takeData[i] <- sum(binData)
	}
	takeData
}