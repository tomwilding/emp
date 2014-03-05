mySdDiff <- function(arr, mean) {
	diffArr <- diff(arr,1)
	mySd(diffArr, mean)
}