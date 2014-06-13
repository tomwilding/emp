mape <- function(preds, data) {
	# Ensure no zeros in data for division
	preds[data == 0] <- 0
	data[data == 0] <- 1
	mape <- sum(abs((preds - data) / data)) / length(data)
} 