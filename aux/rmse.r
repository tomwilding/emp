rmse <- function(preds, data) {
	sqrt(ssError(preds, data) / length(data))
}