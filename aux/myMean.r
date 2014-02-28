myMean <- function(arr) {
	total <- 0
	for (i in 1:length(arr)) {
		total <- total + arr[i]
	}
	total / length(arr)
}