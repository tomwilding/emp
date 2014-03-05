myMeanDiff <- function(arr) {
	diffArr <- diff(arr,1)
	myMean(abs(diffArr))
}