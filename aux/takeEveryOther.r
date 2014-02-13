takeEveryOther <- function(data) {
	takeData <- c()
 	for (i in seq(from=1, to=length(data), by=1)) {
 		if ((i %% 2) == 0) {
 			takeData <- cbind(takeData, data[i])
 		}
	}
	takeData
}