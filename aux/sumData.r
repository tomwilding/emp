sumData <- function(data, n) {
	takeData <- c()
	for (i in 1:(length(data)%/%n)) {
		si <- ((((i-1)*n)+1))
		ei <- (i*n)
		binData <- data[si:ei]
		takeData[i] <- sum(binData)
	}
	takeData
}