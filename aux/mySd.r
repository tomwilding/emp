mySd <- function(arr, m) {
	total <- 0
	for (i in 1:length(arr)) {
		total <- total + (arr[i] - m)^2
	}
	sqrt(total / length(arr))
}