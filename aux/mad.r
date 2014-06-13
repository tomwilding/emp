mad <- function(preds, data) {
	relErrors <- sum(abs(preds - data))
	mad <- relErrors / length(data)
}